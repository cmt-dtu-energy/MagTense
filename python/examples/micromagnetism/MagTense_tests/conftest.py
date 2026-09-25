# The runnable examples in this directory are named <name>_test.py, which matches
# pytest's default "*_test.py" collection pattern. They are scripts driven by
# testMagTenseFunctions.py rather than pytest modules, and importing them at
# collection time pulls in magtense for no benefit, so keep them out. The real
# pytest suites here are all named test_*.py and are unaffected.
collect_ignore_glob = ["*_test.py"]
collect_ignore = ["fmm_test/run_fmm_test.py"]
