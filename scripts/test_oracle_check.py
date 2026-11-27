"""Regression tests for strict cross-platform reference fixture comparison."""
import unittest
from generate_oracle import fixture_equivalent


class OracleComparisonTest(unittest.TestCase):
    baseline = "close(c.log_density(0.3, 0.7), -0.067874571858002, 2.0e-8)"

    def test_linux_libm_last_digit(self):
        self.assertTrue(fixture_equivalent(self.baseline,
            self.baseline.replace("858002", "858001")))

    def test_formatting_and_trailing_comma(self):
        self.assertTrue(fixture_equivalent(self.baseline,
            "close(\n c.log_density(0.3, 0.7,), -0.067874571858002, 2.0e-8,\n)"))

    def test_reference_drift_rejected(self):
        self.assertFalse(fixture_equivalent(self.baseline,
            self.baseline.replace("-0.067874571858002", "-0.067874561858002")))

    def test_input_drift_rejected(self):
        self.assertFalse(fixture_equivalent(self.baseline,
            self.baseline.replace("0.3", "0.300000000000001")))

    def test_tolerance_drift_rejected(self):
        self.assertFalse(fixture_equivalent(self.baseline,
            self.baseline.replace("2.0e-8", "2.000000000000001e-8")))

    def test_function_drift_rejected(self):
        self.assertFalse(fixture_equivalent(self.baseline,
            self.baseline.replace("log_density", "density")))

    def test_data_or_missing_assertion_rejected(self):
        self.assertFalse(fixture_equivalent("let x = [0.3]\n" + self.baseline,
            "let x = [0.300000000000001]\n" + self.baseline))
        self.assertFalse(fixture_equivalent(self.baseline, ""))
        self.assertFalse(fixture_equivalent("", ""))


if __name__ == "__main__":
    unittest.main()
