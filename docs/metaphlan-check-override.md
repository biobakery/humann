# Experimental MetaPhlAn check override

`--override-metaphlan-check` disables MetaPhlAn version checks for the command in which it is specified. It defaults to off. Enabling it emits a `CRITICAL` experimental-use warning explaining that users must understand the compatibility risks and that HUMAnN developers take no responsibility for failures or incorrect results caused by the override.

## Runtime

```sh
humann --input reads.fastq --output results \
    --taxonomic-profile precomputed_profile.tsv \
    --override-metaphlan-check
```

At runtime, the flag skips the `metaphlan --version` database-version probe before prescreening. The warning is written to the HUMAnN run log and stderr even without `--verbose`, and the setting is included in logged run configuration. The flag must be supplied for each run.

When `--taxonomic-profile` is supplied, HUMAnN uses that file instead of executing MetaPhlAn and does not require a MetaPhlAn executable. This also works without the override for a compatible precomputed profile. It supports containers that contain HUMAnN and its alignment dependencies but omit MetaPhlAn.

The override does not remove the requirement for MetaPhlAn when HUMAnN will actually run prescreening. Profile readability, database-header checks, required columns, and other pipeline dependency checks remain enabled. It does not make an incompatible profile format readable.

## Source installation

The function `check_metaphlan_version()` belongs to the custom source install command, where it normally runs before and after installation. The source installer accepts the same flag:

```sh
python setup.py install --override-metaphlan-check
```

With this option, neither call to `check_metaphlan_version()` runs, and a `CRITICAL` experimental-use warning is emitted. Other installation steps still run. This is an option of the existing custom `setup.py install` command; it is not a top-level `pip install` option and does not itself control dependency installation.

Install-time and runtime overrides are independent. An install-time override does not silently disable checks in subsequent HUMAnN runs.

## Focused validation

The focused regression suite verifies opt-in parsing, CRITICAL logging, runtime configuration reset, default version probing, skipped probing under override, unchanged actual execution, precomputed profiles without MetaPhlAn, required executable checks when prescreening is needed, preserved profile header validation, and both source-install version checks.

```sh
python -m unittest humann.tests.basic_tests_metaphlan_override
```

These checks use mocked external commands; they do not run alignments or install dependencies. Full container and end-to-end pipeline testing remains deferred.
