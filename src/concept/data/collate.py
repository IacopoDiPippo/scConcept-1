import os
import numpy as np
import pandas as pd
from typing import Dict, List
from torch.utils.data import Dataset, default_collate, default_convert, get_worker_info
import torch.distributed as dist
from lamin_dataloader import BaseCollate
from abc import ABC, abstractmethod
import random
from pathlib import Path
import re


class Collate(BaseCollate):
    def __init__(self, 
                 tokenizer, 
                 panels_path,
                 max_tokens, 
                 min_tokens=0, 
                 split_input=True, 
                 variable_size=False, 
                 gene_sampling_strategy='top-nonzero',
                 panel_selection='random', 
                 panel_selection_mixed_prob=1.0, 
                 panel_filter_regex='.*', 
                 panel_size_min=None, 
                 panel_size_max=None,
                 panel_overlap=False,
                 panel_max_drop_rate=None, 
                 feature_max_drop_rate=None,
                 model_speed_sanity_check=False,
                 probabilistic_panel_sampling: bool = False,
                 probabilistic_panel_csv: str = None,
                 finetune_panels: bool = False,
                 finetune_selection_mixed_prob: float = 0.25,
                 superselected_panel_1: str = None,
                 superselected_panel_2: str = None,
                 ):
        super().__init__(PAD_TOKEN=tokenizer.PAD_TOKEN, 
                         max_tokens=max_tokens, 
                         gene_sampling_strategy=gene_sampling_strategy, 
        )
        
        self.tokenizer = tokenizer
        self.split_input = split_input
        self.panel_selection = panel_selection
        self.panel_selection_mixed_prob = panel_selection_mixed_prob
        self.finetune_panels = finetune_panels
        self.finetune_selection_mixed_prob = finetune_selection_mixed_prob

        # Load panels for any mode that needs them
        if self.panel_selection not in ('random',) or finetune_panels:
            self.panels_dir = Path(panels_path)
            self.panel_names = [panel_name for panel_name in os.listdir(self.panels_dir) if re.search(panel_filter_regex, panel_name) and panel_name.endswith('.csv')]
            print("Available panels:", self.panel_names)
            self.panels = [self.tokenizer.encode(pd.read_csv(self.panels_dir / panel_name)['Ensembl_ID'].values) 
                       for panel_name in self.panel_names]
            for i in range(len(self.panels)):
                print(f'Panel {self.panel_names[i]} size: {len(self.panels[i])} genes')

        # --- superselected mode: two explicitly named panels ---
        self.superselected_panel_1 = superselected_panel_1
        self.superselected_panel_2 = superselected_panel_2

        if self.panel_selection == 'superselected':
            assert superselected_panel_1 is not None and superselected_panel_2 is not None, \
                "panel_selection='superselected' requires both superselected_panel_1 and superselected_panel_2 to be set"

            self._superselected_idx_1 = self._resolve_superselected_panel(superselected_panel_1, panel_num=1)
            print(f"Superselected panel 1: '{self.panel_names[self._superselected_idx_1]}' ({len(self.panels[self._superselected_idx_1])} genes)")

            self._superselected_idx_2 = self._resolve_superselected_panel(superselected_panel_2, panel_num=2)
            print(f"Superselected panel 2: '{self.panel_names[self._superselected_idx_2]}' ({len(self.panels[self._superselected_idx_2])} genes)")

        self.panel_size_min = panel_size_min
        self.panel_size_max = panel_size_max
        self.panel_overlap = panel_overlap
        self.panel_max_drop_rate = panel_max_drop_rate
        self.gene_sampling_strategy = gene_sampling_strategy
        assert self.gene_sampling_strategy in ['random', 'top', 'random-nonzero', 'top-nonzero'], 'gene_sampling_strategy must be one of "random", "top", "random-nonzero", "top-nonzero"'
        self.feature_max_drop_rate = feature_max_drop_rate
        self.device_num = dist.get_rank() if dist.is_initialized() else 0
        self._rng = None

        self.probabilistic_panel_sampling = probabilistic_panel_sampling

        if self.probabilistic_panel_sampling:
            assert probabilistic_panel_csv is not None, \
                "probabilistic_panel_csv must be provided if probabilistic_panel_sampling=True"

            df = pd.read_csv(probabilistic_panel_csv)
            print("Probabilistic weighted approach for panel sampling enabled.")
            # encode Ensembl IDs → tokens
            tokens = self.tokenizer.encode(df["Ensembl_ID"].values)

            probs = df["score"].values.astype(np.float64)

            self.probabilistic_tokens = np.array(tokens, dtype=np.int64)
            self.probabilistic_probs = np.array(probs, dtype=np.float64)

    def _sample_probabilistic_panel(self, available_tokens, panel_size):
        """
        Sampling without replacement using predefined probabilities.
        Only tokens present in available_tokens are considered.
        """

        mask = np.isin(self.probabilistic_tokens, available_tokens)

        tokens = self.probabilistic_tokens[mask]
        probs = self.probabilistic_probs[mask]

        if len(tokens) == 0:
            # fallback: uniform random
            return self.rng.choice(available_tokens, panel_size, replace=False)

        probs = probs / probs.sum()

        panel = self.rng.choice(
            tokens,
            size=min(panel_size, len(tokens)),
            replace=False,
            p=probs
        )
        return panel


    # This is crucial when running multiple GPUs. 
    # It ensures that the random number generator is the same for each worker.
    @property
    def rng(self):
        if self._rng is None:
            if self.split_input:
                worker_info = get_worker_info()
                if worker_info: # In case of multi-process data loading
                    print(f'Device: {self.device_num}, Worker {worker_info.id} / {worker_info.num_workers}, seed: {worker_info.seed}')
                    self._rng = np.random.default_rng(seed=42 + worker_info.id)
                else:
                    self._rng = np.random.default_rng(42)    
            else:
                self._rng = np.random.default_rng(42)
        return self._rng


    def shared_feature_stats(self, batch):
        num_shared_featrues = []
        for i in range(len(batch)):
            for j in range(len(batch)):
                if i < j :
                    num_shared_featrues.append(np.intersect1d(batch[i]['tokens'], batch[j]['tokens']).size / np.union1d(batch[i]['tokens'], batch[j]['tokens']).size)
        
        print(f'Average % of shared features: %{np.median(num_shared_featrues)*100:.3f}')

    def log_int_samping(self, low, high):
        if low == high:
            return low
        randint = int(np.exp2(self.rng.uniform(np.log2(low), np.log2(high))))
        return max(min(randint, high), low)
    
    def int_samping(self, low, high):
        if low == high:
            return low
        randint = int(self.rng.uniform(low, high))
        return max(min(randint, high), low)
    
    
    def adapt_batch_size(self, batch,  batch_1, batch_2):
        max_length_1 = max(len(item['tokens']) for item in batch_1)
        max_length_2 = max(len(item['tokens']) for item in batch_2)
        
        if max_length_1 > 1000 or max_length_2 > 1000:
            new_batch_size = max(1, len(batch_1) // 4)
            batch = batch[:new_batch_size]
            batch_1 = batch_1[:new_batch_size]
            batch_2 = batch_2[:new_batch_size]
        
        return batch, batch_1, batch_2

    def _resolve_superselected_panel(self, name: str, panel_num: int) -> int:
        """
        Resolve a superselected panel name to its index in self.panel_names.

        The name must match the panel filename (without the .csv extension) exactly.
        If no exact match is found, raises a ValueError listing all available panels
        so the user can correct the config immediately.
        """
        # Strip .csv suffix from the provided name if the user accidentally included it
        lookup = name.removesuffix('.csv')

        # Try exact match against bare names (filename without extension)
        bare_names = [n.removesuffix('.csv') for n in self.panel_names]
        if lookup in bare_names:
            return bare_names.index(lookup)

        # No match — build a helpful error
        available = '\n  '.join(sorted(bare_names))
        raise ValueError(
            f"\n\n[superselected_panel_{panel_num}] Panel '{name}' was not found in the panels directory.\n"
            f"Available panels (filename without .csv):\n  {available}\n\n"
            f"Please update your config to use one of the names above exactly."
        )

    def _get_predesigned_panel(self, batch):

        i = self.rng.integers(0, len(self.panels))
        panel = self.panels[i]
        if self.panel_max_drop_rate is not None and self.panel_max_drop_rate > 0:
            panel_max_drop_rate = self.rng.uniform(0, self.panel_max_drop_rate)
            drop_mask = self.rng.uniform(size=len(panel)) > panel_max_drop_rate
            panel = panel[drop_mask]
        return panel, self.panel_names[i]

    def _get_panel_by_index(self, batch, panel_idx):
        """Return a specific panel by its index (used for superselected mode)."""
        panel = self.panels[panel_idx]
        if self.panel_max_drop_rate is not None and self.panel_max_drop_rate > 0:
            panel_max_drop_rate = self.rng.uniform(0, self.panel_max_drop_rate)
            drop_mask = self.rng.uniform(size=len(panel)) > panel_max_drop_rate
            panel = panel[drop_mask]
        return panel, self.panel_names[panel_idx]
    
    def __call__(self, batch):
        n_tokens = len(batch[0]['tokens'])
        permute = self.rng.permutation(n_tokens)
        batch_permute = [{'tokens': item['tokens'][permute], 
                          'values': item['values'][permute], 
                          } for item in batch]
        
        if self.split_input:
            
            n_tokens = len(batch_permute[0]['tokens'])
            panel_indices = np.arange(n_tokens)

            panel_name_1, panel_name_2 = 'random', 'random'
            panel_overlap = self.rng.uniform() <= float(self.panel_overlap)

            # -------------------------------------------------------
            # SUPERSELECTED: both panels are explicitly named
            # -------------------------------------------------------
            if self.panel_selection == 'superselected':
                panel_1_tokens, panel_name_1 = self._get_panel_by_index(batch_permute, self._superselected_idx_1)
                panel_idx_1 = np.where(np.isin(batch_permute[0]['tokens'], panel_1_tokens))[0]

                panel_2_tokens, panel_name_2 = self._get_panel_by_index(batch_permute, self._superselected_idx_2)
                panel_idx_2 = np.where(np.isin(batch_permute[0]['tokens'], panel_2_tokens))[0]
                # Note: superselected panels CAN overlap (they are fixed biology panels)

            # -------------------------------------------------------
            # RANDOM / MIXED / PRESELECTED (existing logic)
            # -------------------------------------------------------
            elif self.panel_selection == 'random' or (self.panel_selection == 'mixed' and self.rng.uniform() <= self.panel_selection_mixed_prob) or n_tokens < 10_000:
                n_tokens_available = n_tokens if panel_overlap else max((n_tokens - self.panel_size_min), 0)
                panel_size_1 = self.log_int_samping(min(self.panel_size_min, n_tokens_available), min(self.panel_size_max, n_tokens_available))
                if self.probabilistic_panel_sampling:
                    panel_tokens = self._sample_probabilistic_panel(
                        batch_permute[0]["tokens"],
                        panel_size_1
                    )
                    panel_idx_1 = np.where(
                        np.isin(batch_permute[0]["tokens"], panel_tokens)
                    )[0]
                else:
                    panel_idx_1 = self.rng.choice(
                        panel_indices,
                        panel_size_1,
                        replace=False
                    )

            else:
                panel, panel_name_1 = self._get_predesigned_panel(batch_permute)
                panel_idx_1 = np.where(np.isin(batch_permute[0]['tokens'], panel))[0]
                panel_size_1 = len(panel_idx_1)

            # -------------------------------------------------------
            # Panel 2 selection (only for non-superselected modes)
            # -------------------------------------------------------
            if self.panel_selection != 'superselected':
                if self.finetune_panels and self.rng.uniform() <= self.finetune_selection_mixed_prob:
                    panel, panel_name_2 = self._get_predesigned_panel(batch_permute)
                    panel_idx_2 = np.where(np.isin(batch_permute[0]['tokens'], panel))[0]
                    panel_size_2 = len(panel_idx_2)

                elif panel_overlap: 
                    panel_size_2 = self.log_int_samping(
                        min(self.panel_size_min, n_tokens),
                        min(self.panel_size_max, n_tokens)
                    )

                    if self.probabilistic_panel_sampling:
                        panel_tokens_2 = self._sample_probabilistic_panel(
                            batch_permute[0]["tokens"],
                            panel_size_2
                        )
                        panel_idx_2 = np.where(
                            np.isin(batch_permute[0]["tokens"], panel_tokens_2)
                        )[0]
                    else:
                        panel_idx_2 = self.rng.choice(
                            panel_indices,
                            panel_size_2,
                            replace=False
                        )

                else:
                    panel_size_2 = self.log_int_samping(
                        min(self.panel_size_min, n_tokens - panel_size_1),
                        min(self.panel_size_max, n_tokens - panel_size_1)
                    )

                    available_tokens = batch_permute[0]["tokens"][~np.isin(
                        batch_permute[0]["tokens"],
                        batch_permute[0]["tokens"][panel_idx_1]
                    )]

                    if self.probabilistic_panel_sampling:
                        panel_tokens_2 = self._sample_probabilistic_panel(
                            available_tokens,
                            panel_size_2
                        )
                        panel_idx_2 = np.where(
                            np.isin(batch_permute[0]["tokens"], panel_tokens_2)
                        )[0]
                    else:
                        panel_idx_2 = self.rng.choice(
                            np.setdiff1d(panel_indices, panel_idx_1, assume_unique=True),
                            panel_size_2,
                            replace=False
                        )

                    assert np.intersect1d(panel_idx_1, panel_idx_2).size == 0, 'Panels overlap'
            
                        
            batch_1 = [{'tokens': item['tokens'][panel_idx_1], 'values': item['values'][panel_idx_1]} for item in batch_permute]
            batch_2 = [{'tokens': item['tokens'][panel_idx_2], 'values': item['values'][panel_idx_2]} for item in batch_permute]

            # The following can be optimized by only passing one panel per batch
            panel_1 = [item['tokens'] for item in batch_1]
            panel_2 = [item['tokens'] for item in batch_2]

            batch_1 = [self.select_features(item, self.feature_max_drop_rate) for item in batch_1]
            batch_2 = [self.select_features(item, self.feature_max_drop_rate) for item in batch_2]

            max_lenght_1 = max([len(item['tokens']) for item in batch_1])
            max_lenght_1 = min(max_lenght_1, self.max_tokens - 1) # todo
            max_lenght_2 = max([len(item['tokens']) for item in batch_2])
            max_lenght_2 = min(max_lenght_2, self.max_tokens - 1) # todo
            
            
            batch_1 = [self.resize_and_pad(item, max_lenght_1) for item in batch_1]
            batch_2 = [self.resize_and_pad(item, max_lenght_2) for item in batch_2]
            
            tokens_1 = [item['tokens'].astype(np.int64) for item in batch_1]
            values_1 = [item['values'].astype(np.float32) for item in batch_1]
            tokens_2 = [item['tokens'].astype(np.int64) for item in batch_2]
            values_2 = [item['values'].astype(np.float32) for item in batch_2]
            
            return {'tokens_1': default_collate(tokens_1),
                    'values_1': default_collate(values_1),
                    'tokens_2': default_collate(tokens_2),
                    'values_2': default_collate(values_2),
                    'panel_1' : default_collate(panel_1),
                    'panel_2' : default_collate(panel_2),
                    'panel_name_1': panel_name_1,
                    'panel_name_2': panel_name_2,
                    **{key: default_collate([item[key] for item in batch]) for key in batch[0].keys() if key not in ['tokens', 'values']}
            }

        else:
            batch_ = [self.select_features(item) for item in batch_permute]
            
            max_lenght = max([len(item['tokens']) for item in batch_])
            
            max_lenght = min(max_lenght, self.max_tokens - 1)

            batch_ = [self.resize_and_pad(item, max_lenght) for item in batch_]
            
            tokens, values = [item['tokens'].astype(np.int64) for item in batch_], [item['values'].astype(np.float32) for item in batch_]
        
            return {'tokens': default_collate(tokens),
                    'values': default_collate(values),
                    **{key: default_collate([item[key] for item in batch]) for key in batch[0].keys() if key not in ['tokens', 'values']}
            }