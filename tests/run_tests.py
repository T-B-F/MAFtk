#!/usr/bin/env python
"""PyUnit-based tests for MAFtk.
"""

from __future__ import print_function

# standard modules
import sys
import os
import re
import getopt
import time
import traceback
import unittest
import distutils.util
import gc
from io import BytesIO

try:
    from StringIO import StringIO  # Python 2 (byte strings)
except ImportError:
    from io import StringIO  # Python 3 (unicode strings)

def is_pypy():
    import platform
    try:
        if platform.python_implementation() == 'PyPy':
            return True
    except AttributeError:
        # New in Python 2.6, not in Jython yet either
        pass
    return False

def is_mafio():
    try:
        from Bio.AlignIO import MafIO
        del MafIO
        return True
    except ImportError:
        return False

def is_intervaltree():
    try:
        import intervaltree
        del intervaltree
        return True
    except ImportError:
        return False

# The default verbosity (not verbose)
VERBOSITY = 0

system_lang = os.environ.get('LANG', 'C')  # Cache this

def main(argv):
    """Run tests, return number of failures (integer)."""
    # insert our paths in sys.path:
    # ../build/lib.*
    # ..
    # Q. Why this order?
    # A. To find the C modules (which are in ../build/lib.*/Bio)
    # Q. Then, why ".."?
    # A. Because Martel may not be in ../build/lib.*
    test_path = sys.path[0] or "."
    source_path = os.path.abspath("%s/../lib/" % test_path)
    sys.path.insert(1, source_path)
    build_path = os.path.abspath("%s/../build/lib.%s-%s" % (
        test_path, distutils.util.get_platform(), sys.version[:3]))
    if os.access(build_path, os.F_OK):
        sys.path.insert(1, build_path)

    # Using "export LANG=C" (which should work on Linux and similar) can
    # avoid problems detecting optional command line tools on
    # non-English OS (we may want 'command not found' in English).
    # HOWEVER, we do not want to change the default encoding which is
    # rather important on Python 3 with unicode.
    # lang = os.environ['LANG']

    # get the command line options
    try:
        opts, args = getopt.getopt(argv, 'gv', ["generate", "verbose", "help"])
    except getopt.error as msg:
        print(msg)
        print(__doc__)
        return 2

    verbosity = VERBOSITY

    # deal with the options
    for o, a in opts:
        if o == "--help":
            print(__doc__)
            return 0
        if o == "-g" or o == "--generate":
            if len(args) > 1:
                print("Only one argument (the test name) needed for generate")
                print(__doc__)
                return 2
            elif len(args) == 0:
                print("No test name specified to generate output for.")
                print(__doc__)
                return 2
            # strip off .py if it was included
            if args[0][-3:] == ".py":
                args[0] = args[0][:-3]

            test = ComparisonTestCase(args[0])
            test.generate_output()
            return 0

        if o == "-v" or o == "--verbose":
            verbosity = 2

    # deal with the arguments, which should be names of tests to run
    for arg_num in range(len(args)):
        # strip off the .py if it was included
        if args[arg_num][-3:] == ".py":
            args[arg_num] = args[arg_num][:-3]

    print("Python version: %s" % sys.version)
    print("Operating system: %s %s" % (os.name, sys.platform))

    # run the tests
    runner = TestRunner(args, verbosity)
    return runner.run()


class ComparisonTestCase(unittest.TestCase):
    """Run a print-and-compare test and check it against expected output."""

    def __init__(self, name, output=None):
        """Initialize with the test to run.

        Arguments:
        o name - The name of the test. The expected output should be
          stored in the file output/name.
        o output - The output that was generated when this test was run.
        """
        unittest.TestCase.__init__(self)
        self.name = name
        self.output = output

    def shortDescription(self):
        return self.name

    def runTest(self):
        # check the expected output to be consistent with what
        # we generated
        outputdir = os.path.join(TestRunner.testdir, "output")
        outputfile = os.path.join(outputdir, self.name)
        try:
            if sys.version_info[0] >= 3:
                # Python 3 problem: Can't use utf8 on output/test_geo
                # due to micro (\xb5) and degrees (\xb0) symbols
                # Also universal new lines mode deprecated on Python 3
                expected = open(outputfile, encoding="latin")
            else:
                expected = open(outputfile, "rU")
        except IOError:
            self.fail("Warning: Can't open %s for test %s" %
                      (outputfile, self.name))

        self.output.seek(0)
        # first check that we are dealing with the right output
        # the first line of the output file is the test name
        expected_test = expected.readline().strip()

        if expected_test != self.name:
            expected.close()
            raise ValueError("\nOutput:   %s\nExpected: %s"
                             % (self.name, expected_test))

        # now loop through the output and compare it to the expected file
        while True:
            expected_line = expected.readline()
            output_line = self.output.readline()

            # stop looping if either of the info handles reach the end
            if not(expected_line) or not(output_line):
                # make sure both have no information left
                assert expected_line == '', "Unread: %s" % expected_line
                assert output_line == '', "Extra output: %s" % output_line
                break

            # normalize the newlines in the two lines
            expected_line = expected_line.strip("\r\n")
            output_line = output_line.strip("\r\n")

            # if the line is a PyUnit time output like:
            # Ran 2 tests in 0.285s
            # ignore it, so we don't have problems with different running times
            if re.compile("^Ran [0-9]+ tests? in ").match(expected_line):
                pass
            # otherwise make sure the two lines are the same
            elif expected_line != output_line:
                expected.close()
                raise ValueError("\nOutput  : %s\nExpected: %s"
                                 % (repr(output_line), repr(expected_line)))
        expected.close()

    def generate_output(self):
        """Generate the golden output for the specified test.
        """
        outputdir = os.path.join(TestRunner.testdir, "output")
        outputfile = os.path.join(outputdir, self.name)

        output_handle = open(outputfile, 'w')

        # write the test name as the first line of the output
        output_handle.write(self.name + "\n")

        # remember standard out so we can reset it after we are done
        save_stdout = sys.stdout
        try:
            # write the output from the test into a string
            sys.stdout = output_handle
            __import__(self.name)
        finally:
            output_handle.close()
            # return standard out to its normal setting
            sys.stdout = save_stdout


class TestRunner(unittest.TextTestRunner):

    if __name__ == '__main__':
        file = sys.argv[0]
    else:
        file = __file__
    testdir = os.path.abspath(os.path.dirname(file) or os.curdir)

    def __init__(self, tests=(), verbosity=0):
        # if no tests were specified to run, we run them all
        self.tests = tests
        if not self.tests:
            # Make a list of all applicable test modules.
            names = os.listdir(TestRunner.testdir)
            for name in names:
                if name[:5] == "test_" and name[-3:] == ".py":
                    self.tests.append(name[:-3])
            self.tests.sort()
        stream = StringIO()
        unittest.TextTestRunner.__init__(self, stream,
                                         verbosity=verbosity)

    def runTest(self, name):
        result = self._makeResult()
        output = StringIO()
        # Restore the language and thus default encoding (in case a prior
        # test changed this, e.g. to help with detecting command line tools)
        global system_lang
        os.environ['LANG'] = system_lang
        # Always run tests from the Tests/ folder where run_tests.py
        # should be located (as we assume this with relative paths etc)
        os.chdir(self.testdir)
        try:
            stdout = sys.stdout
            sys.stdout = output
            if name.startswith("test_"):
                sys.stderr.write("%s ... " % name)
                # It's either a unittest or a print-and-compare test
                loader = unittest.TestLoader()
                suite = loader.loadTestsFromName(name)
                if hasattr(loader, "errors") and loader.errors:                       
                    # Looks like a real failure
                    sys.stderr.write("loading tests failed:\n")
                    for msg in loader.errors:
                        sys.stderr.write("%s\n" % msg)
                    return False
                if suite.countTestCases() == 0:
                    # This is a print-and-compare test instead of a
                    # unittest-type test.
                    test = ComparisonTestCase(name, output)
                    suite = unittest.TestSuite([test])
            suite.run(result)
            if self.testdir != os.path.abspath("."):
                sys.stderr.write("FAIL\n")
                result.stream.write(result.separator1 + "\n")
                result.stream.write("ERROR: %s\n" % name)
                result.stream.write(result.separator2 + "\n")
                result.stream.write("Current directory changed\n")
                result.stream.write("Was: %s\n" % self.testdir)
                result.stream.write("Now: %s\n" % os.path.abspath("."))
                os.chdir(self.testdir)
                if not result.wasSuccessful():
                    result.printErrors()
                return False
            elif result.wasSuccessful():
                sys.stderr.write("ok\n")
                return True
            else:
                sys.stderr.write("FAIL\n")
                result.printErrors()
            return False
        except Exception as msg:
            # This happened during the import
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write("ERROR: %s\n" % name)
            result.stream.write(result.separator2 + "\n")
            result.stream.write(traceback.format_exc())
            return False
        except KeyboardInterrupt as err:
            # Want to allow this, and abort the test
            # (see below for special case)
            raise err
        except:
            # This happens in Jython with java.lang.ClassFormatError:
            # Invalid method Code length ...
            sys.stderr.write("ERROR\n")
            result.stream.write(result.separator1 + "\n")
            result.stream.write("ERROR: %s\n" % name)
            result.stream.write(result.separator2 + "\n")
            result.stream.write(traceback.format_exc())
            return False
        finally:
            sys.stdout = stdout
            # Running under PyPy we were leaking file handles...
            gc.collect()

    def run(self):
        """Run tests, return number of failures (integer)."""
        failures = 0
        startTime = time.time()
        for test in self.tests:
            ok = self.runTest(test)
            if not ok:
                failures += 1
        total = len(self.tests)
        stopTime = time.time()
        timeTaken = stopTime - startTime
        sys.stderr.write(self.stream.getvalue())
        sys.stderr.write('-' * 70 + "\n")
        sys.stderr.write("Ran %d test%s in %.3f seconds\n" %
                         (total, total != 1 and "s" or "", timeTaken))
        sys.stderr.write("\n")
        if failures:
            sys.stderr.write("FAILED (failures = %d)\n" % failures)
        return failures


if __name__ == "__main__":
    errors = main(sys.argv[1:])
    if errors:
        # Doing a sys.exit(...) isn't nice if run from IDLE...
        sys.exit(1)
