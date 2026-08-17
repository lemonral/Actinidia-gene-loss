# Tests

Run the core gene-loss smoke test without any external scientific package:

```bash
PYTHONPATH=src python3 -m unittest discover -s tests -v
```

The test creates a tiny genome/GFF, normalizes a 13-column SynOrths file,
calls a bracketed candidate, parses a 12-column negative-strand tBLASTX hit,
creates the classification and runs the five-bin spatial summary. It also
checks the explicit legacy six-column parser path.

It does not test SynOrths, BLAST+, large-genome memory use, SciPy statistics, or
Matplotlib rendering. Run those integration tests on a small representative
dataset before regenerating study outputs.
