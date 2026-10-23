# Limitations and release notes

This repository implements a source-list-to-database pipeline; it does not distribute the original complete OCVRC master-catalog research workspace. Survey completeness depends on sky coverage, processing flags and service availability. A positive fitted flux alone is not a validated astronomical detection.

The 2026-10-23 six-source test completed all pipeline stages and produced a validated 82-row measurement table. SIMBAD source-ingestion outputs were **reused** from an earlier successful run; thus that test did not exercise a fresh live SIMBAD query or an initial empty-cache build.

Before public release, verify a fresh-input run, package dependencies on a clean environment and determine the project's licensing terms. The publicly hosted older MkDocs site may contain additional historical explanatory text not supplied in this archive.
