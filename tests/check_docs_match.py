#!/usr/bin/env python3
"""Check that docs/source/test_suite.rst lists exactly the tests that exist.

The page drifted from the suite the moment a test was added, and nothing
noticed.  Comparing the two mechanically is cheap, and a table that claims to
list the tests is worth no more than its accuracy.
"""
import io, os, re, subprocess, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
PAGE = os.path.join(ROOT, 'docs', 'source', 'test_suite.rst')


def tests_in_suite():
    out = subprocess.run([sys.executable, os.path.join(HERE, 'run_tests.py'), '--list'],
                         capture_output=True, text=True, check=True).stdout
    return {line.split()[0] for line in out.splitlines() if line.strip()}


def tests_on_page():
    text = io.open(PAGE, encoding='utf-8').read()
    # rows of the table read "   * - ``name``"
    return set(re.findall(r'^\s*\* - ``(\w+)``\s*$', text, re.M))


def main():
    have, listed = tests_in_suite(), tests_on_page()
    missing = sorted(have - listed)
    extra = sorted(listed - have)
    if not missing and not extra:
        print('test_suite.rst lists all %d tests' % len(have))
        return 0
    if missing:
        print('in the suite but not on the page: %s' % ', '.join(missing))
    if extra:
        print('on the page but not in the suite: %s' % ', '.join(extra))
    print('\nEdit docs/source/test_suite.rst so the table matches.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
