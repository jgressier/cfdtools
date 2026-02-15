import functools
import tracemalloc
import time
from cfdtools import log

def lazyprop(fn):
    attr_name = '_cache_' + fn.__name__

    @property
    @functools.wraps(fn)
    def _lazyprop(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fn(self))
        return getattr(self, attr_name)

    return _lazyprop


def cpu_time_decorator(func):
    """decorator to print the elapsed time of a function
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        start_time = time.perf_counter()  # Start the timer
        result = func(*args, **kwargs)    # Execute the function
        end_time = time.perf_counter()    # End the timer
        elapsed_time = end_time - start_time  # Calculate elapsed time

        # Print the function name and elapsed time with right alignment
        log.debug(f"{func.__name__:>75}: {elapsed_time:.1f}s")

        return result
    return wrapper

class TraceMemory:
    def __init__(self):
        self._snapshot = None

    def __enter__(self):
        tracemalloc.start()
        self._snapshot = tracemalloc.take_snapshot()

    def __exit__(self, *exitoptions):
        current, peak = tracemalloc.get_traced_memory()
        # snapshot = tracemalloc.take_snapshot()
        # top_stats = snapshot.compare_to(self._snapshot, 'lineno')
        # print("[ Top 10 ]")
        # for stat in top_stats[:10]:
        #     print(stat)
        log.debug(f"  Current Memory Usage: {current / 1e6:.1f} MB")
        log.debug(f"     Peak Memory Usage: {peak / 1e6:.1f} MB")
        # Stop tracing memory allocation
        tracemalloc.stop()
