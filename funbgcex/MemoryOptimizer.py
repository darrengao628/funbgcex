#!/usr/bin/env python3

"""
Memory optimization utilities for FunBGCeX
This module provides utilities for memory-efficient processing of large datasets.
"""

import os
import psutil
import hashlib
import pickle
import time
from functools import lru_cache
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Generator, Any
import pandas as pd
import dask.dataframe as dd
from dask.diagnostics import ProgressBar
import logging

logger = logging.getLogger(__name__)


class MemoryMonitor:
    """
    Enhanced memory monitor with configurable limits and warnings.
    """
    def __init__(self, max_memory_gb: Optional[float] = None, 
                 warning_threshold: float = 0.8, 
                 check_interval: int = 5):
        """
        Initialize memory monitor.
        
        Args:
            max_memory_gb: Maximum memory usage in GB (None for no limit)
            warning_threshold: Threshold for memory warnings (0.0-1.0)
            check_interval: Interval in seconds between memory checks
        """
        self.max_memory_gb = max_memory_gb
        self.warning_threshold = warning_threshold
        self.check_interval = check_interval
        self.process = psutil.Process()
        self.initial_memory = psutil.virtual_memory().used / (1024**3)  # GB
        self.monitoring = False
        self.last_warning_time = 0
        
    def start_monitoring(self):
        """Start memory monitoring."""
        if self.max_memory_gb is not None:
            self.monitoring = True
            logger.info(f"Memory monitoring started. Limit: {self.max_memory_gb} GB")
            
    def stop_monitoring(self):
        """Stop memory monitoring."""
        self.monitoring = False
        logger.info("Memory monitoring stopped")
        
    def check_memory(self) -> Tuple[float, float]:
        """
        Check current memory usage.
        
        Returns:
            Tuple of (current_usage_gb, percentage_of_max)
        """
        current_memory = psutil.virtual_memory().used / (1024**3)  # GB
        if self.max_memory_gb is not None:
            percentage = (current_memory / self.max_memory_gb) * 100
        else:
            percentage = 0
        return current_memory, percentage
        
    def is_memory_exceeded(self) -> bool:
        """
        Check if memory usage exceeds the limit.
        
        Returns:
            True if memory usage exceeds the limit, False otherwise
        """
        if self.max_memory_gb is None:
            return False
            
        current_memory, percentage = self.check_memory()
        
        # Log warning if approaching limit
        if (percentage > self.warning_threshold * 100 and 
            time.time() - self.last_warning_time > 60):  # Warn at most once per minute
            logger.warning(f"Memory usage high: {percentage:.1f}% of limit ({current_memory:.2f} GB)")
            self.last_warning_time = time.time()
            
        return percentage > 100
        
    def get_memory_info(self) -> Dict[str, float]:
        """
        Get detailed memory information.
        
        Returns:
            Dictionary with memory information
        """
        mem = psutil.virtual_memory()
        return {
            'total': mem.total / (1024**3),  # GB
            'available': mem.available / (1024**3),  # GB
            'used': mem.used / (1024**3),  # GB
            'percent': mem.percent,
            'initial': self.initial_memory,
            'current': mem.used / (1024**3),
            'max_limit': self.max_memory_gb
        }


class StreamingDataFrame:
    """
    Wrapper for memory-efficient DataFrame operations using Dask.
    """
    def __init__(self, chunk_size: int = 10000, max_memory_gb: Optional[float] = None):
        """
        Initialize streaming DataFrame processor.
        
        Args:
            chunk_size: Number of rows to process at a time
            max_memory_gb: Maximum memory usage in GB
        """
        self.chunk_size = chunk_size
        self.memory_monitor = MemoryMonitor(max_memory_gb)
        
    def read_csv_streaming(self, filepath: str, **kwargs) -> dd.DataFrame:
        """
        Read CSV file in streaming mode using Dask.
        
        Args:
            filepath: Path to CSV file
            **kwargs: Additional arguments for dd.read_csv
            
        Returns:
            Dask DataFrame
        """
        logger.info(f"Reading CSV in streaming mode: {filepath}")
        return dd.read_csv(filepath, blocksize=self.chunk_size * 1024, **kwargs)
        
    def process_dataframe_streaming(self, df: pd.DataFrame, 
                                  processing_func, 
                                  **kwargs) -> pd.DataFrame:
        """
        Process DataFrame in chunks to reduce memory usage.
        
        Args:
            df: Input DataFrame
            processing_func: Function to apply to each chunk
            **kwargs: Additional arguments for processing_func
            
        Returns:
            Processed DataFrame
        """
        logger.info(f"Processing DataFrame in chunks of {self.chunk_size} rows")
        
        # Convert to Dask DataFrame
        ddf = dd.from_pandas(df, npartitions=max(1, len(df) // self.chunk_size))
        
        # Apply processing function
        result_ddf = ddf.map_partitions(processing_func, meta=df, **kwargs)
        
        # Compute with progress bar
        with ProgressBar():
            result_df = result_ddf.compute()
            
        logger.info(f"Processed {len(result_df)} rows")
        return result_df


class SequenceProcessor:
    """
    Memory-efficient sequence processing using generators.
    """
    def __init__(self, cache_size: int = 1000):
        """
        Initialize sequence processor.
        
        Args:
            cache_size: Number of sequences to cache in memory
        """
        self.cache_size = cache_size
        
    @lru_cache(maxsize=1000)
    def get_sequence_hash(self, sequence: str) -> str:
        """
        Get hash of a sequence for caching purposes.
        
        Args:
            sequence: Protein sequence
            
        Returns:
            SHA256 hash of the sequence
        """
        return hashlib.sha256(sequence.encode()).hexdigest()
        
    def read_sequences_generator(self, filepath: str) -> Generator[Tuple[str, str], None, None]:
        """
        Read sequences from FASTA file using a generator.
        
        Args:
            filepath: Path to FASTA file
            
        Yields:
            Tuple of (sequence_id, sequence)
        """
        current_id = None
        current_seq = []
        
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_id is not None:
                        yield current_id, ''.join(current_seq)
                    current_id = line[1:]  # Remove '>'
                    current_seq = []
                else:
                    current_seq.append(line)
                    
        # Yield last sequence
        if current_id is not None:
            yield current_id, ''.join(current_seq)
            
    def process_sequences_streaming(self, filepath: str, 
                                   processing_func,
                                   **kwargs) -> List[Tuple[str, Any]]:
        """
        Process sequences from FASTA file in streaming mode.
        
        Args:
            filepath: Path to FASTA file
            processing_func: Function to apply to each sequence
            **kwargs: Additional arguments for processing_func
            
        Returns:
            List of (sequence_id, result) tuples
        """
        logger.info(f"Processing sequences in streaming mode: {filepath}")
        results = []
        
        for seq_id, sequence in self.read_sequences_generator(filepath):
            # Check memory before processing
            if hasattr(self, 'memory_monitor') and self.memory_monitor.is_memory_exceeded():
                logger.warning("Memory limit exceeded, stopping processing")
                break
                
            result = processing_func(seq_id, sequence, **kwargs)
            results.append((seq_id, result))
            
        logger.info(f"Processed {len(results)} sequences")
        return results


class DiskCache:
    """
    Disk-based caching system for expensive computations.
    """
    def __init__(self, cache_dir: str = ".cache", max_size_mb: int = 1000):
        """
        Initialize disk cache.
        
        Args:
            cache_dir: Directory to store cache files
            max_size_mb: Maximum cache size in MB
        """
        self.cache_dir = Path(cache_dir)
        self.cache_dir.mkdir(exist_ok=True)
        self.max_size_mb = max_size_mb
        self.index_file = self.cache_dir / "index.pkl"
        
        # Load or create cache index
        if self.index_file.exists():
            with open(self.index_file, 'rb') as f:
                self.index = pickle.load(f)
        else:
            self.index = {}
            
    def _get_cache_key(self, func_name: str, args: tuple, kwargs: dict) -> str:
        """
        Generate cache key from function name and arguments.
        
        Args:
            func_name: Name of the function
            args: Function arguments
            kwargs: Function keyword arguments
            
        Returns:
            Cache key string
        """
        # Create a hash of the arguments
        args_str = str(args) + str(sorted(kwargs.items()))
        args_hash = hashlib.md5(args_str.encode()).hexdigest()
        return f"{func_name}_{args_hash}"
        
    def _get_cache_path(self, cache_key: str) -> Path:
        """
        Get file path for cache key.
        
        Args:
            cache_key: Cache key
            
        Returns:
            Path to cache file
        """
        return self.cache_dir / f"{cache_key}.pkl"
        
    def get(self, func_name: str, args: tuple, kwargs: dict) -> Optional[Any]:
        """
        Get cached result.
        
        Args:
            func_name: Name of the function
            args: Function arguments
            kwargs: Function keyword arguments
            
        Returns:
            Cached result or None if not found
        """
        cache_key = self._get_cache_key(func_name, args, kwargs)
        cache_path = self._get_cache_path(cache_key)
        
        if cache_path.exists():
            # Check if cache is still valid
            if cache_key in self.index:
                file_mtime = cache_path.stat().st_mtime
                if file_mtime >= self.index[cache_key]['input_mtime']:
                    with open(cache_path, 'rb') as f:
                        logger.debug(f"Cache hit: {cache_key}")
                        return pickle.load(f)
                else:
                    # Cache is invalid
                    self._remove_cache_entry(cache_key)
                    
        return None
        
    def set(self, func_name: str, args: tuple, kwargs: dict, 
           result: Any, input_mtime: float):
        """
        Set cached result.
        
        Args:
            func_name: Name of the function
            args: Function arguments
            kwargs: Function keyword arguments
            result: Result to cache
            input_mtime: Modification time of input file
        """
        cache_key = self._get_cache_key(func_name, args, kwargs)
        cache_path = self._get_cache_path(cache_key)
        
        # Check cache size and clean up if necessary
        self._cleanup_cache()
        
        # Save result
        with open(cache_path, 'wb') as f:
            pickle.dump(result, f)
            
        # Update index
        self.index[cache_key] = {
            'file_path': str(cache_path),
            'size': cache_path.stat().st_size,
            'input_mtime': input_mtime,
            'access_time': time.time()
        }
        
        # Save index
        with open(self.index_file, 'wb') as f:
            pickle.dump(self.index, f)
            
        logger.debug(f"Cache set: {cache_key}")
        
    def _cleanup_cache(self):
        """Clean up cache if it exceeds size limit."""
        total_size = sum(entry['size'] for entry in self.index.values())
        total_size_mb = total_size / (1024 * 1024)
        
        if total_size_mb > self.max_size_mb:
            # Remove least recently used entries
            sorted_entries = sorted(self.index.items(), 
                                  key=lambda x: x[1]['access_time'])
            
            for cache_key, entry in sorted_entries:
                self._remove_cache_entry(cache_key)
                total_size_mb -= entry['size'] / (1024 * 1024)
                
                if total_size_mb <= self.max_size_mb * 0.8:  # Clean to 80% of limit
                    break
                    
    def _remove_cache_entry(self, cache_key: str):
        """Remove cache entry."""
        if cache_key in self.index:
            cache_path = Path(self.index[cache_key]['file_path'])
            if cache_path.exists():
                cache_path.unlink()
            del self.index[cache_key]
            
    def clear(self):
        """Clear all cache entries."""
        for cache_key in list(self.index.keys()):
            self._remove_cache_entry(cache_key)
        self.index = {}
        with open(self.index_file, 'wb') as f:
            pickle.dump(self.index, f)


def cached_computation(cache_dir: str = ".cache", max_size_mb: int = 1000):
    """
    Decorator for caching function results to disk.
    
    Args:
        cache_dir: Directory to store cache files
        max_size_mb: Maximum cache size in MB
    """
    cache = DiskCache(cache_dir, max_size_mb)
    
    def decorator(func):
        def wrapper(*args, **kwargs):
            # Get input file modification time if available
            input_mtime = 0
            if 'input_file' in kwargs:
                input_path = Path(kwargs['input_file'])
                if input_path.exists():
                    input_mtime = input_path.stat().st_mtime
                    
            # Try to get from cache
            cached_result = cache.get(func.__name__, args, kwargs)
            if cached_result is not None:
                return cached_result
                
            # Compute result
            result = func(*args, **kwargs)
            
            # Cache result
            cache.set(func.__name__, args, kwargs, result, input_mtime)
            
            return result
        return wrapper
    return decorator