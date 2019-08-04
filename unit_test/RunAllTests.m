% run the tests for Travis CI

import matlab.unittest.plugins.CodeCoveragePlugin
import matlab.unittest.plugins.codecoverage.CoberturaFormat
import matlab.unittest.TestRunner

suite = testsuite('IncludeSubfolders', false);
runner = TestRunner.withTextOutput;

% add a coverage report
reportFile = fullfile('..', 'coverage.xml');
reportFormat = CoberturaFormat(reportFile);
plugin = CodeCoveragePlugin.forFolder('..', 'Producing',reportFormat);
runner.addPlugin(plugin);

%% setup the path
addpath ..

%% Run all unit tests in my repository

% Run all unit tests in my repository.
results = runner.run(suite);

% Assert no tests failed.
assert(all(~[results.Failed]));
