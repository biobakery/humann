"""Preserve raw read-count metadata independently of abundance rows."""
import logging
from . import config

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)


def parse_line(line):
    if line.startswith(config.read_count_comment_prefix):
        fields = line[len(config.read_count_comment_prefix):].rstrip("\r\n").split("\t")
        if fields[0] in config.read_count_names:
            return fields[0], fields[1:]
    return None


def comment_lines(read_counts):
    for name in config.read_count_names:
        if name in read_counts:
            values = read_counts[name]
            if not isinstance(values, (list, tuple)):
                values = [values]
            values = [config.read_count_unavailable if value is None else
                      value if isinstance(value, str) else "%.15g" % value for value in values]
            yield config.read_count_comment_prefix + name + "\t" + "\t".join(values)


def log_counts(read_counts):
    for line in comment_lines(read_counts):
        logger.info("%s", line[len(config.read_count_comment_prefix):])


class LegacyProfileWarning:
    """Warn once per ingestion without changing legacy profile values."""

    def __init__(self, path):
        self.path = path
        self.warned = False

    def check(self, line):
        feature = line.split(config.gene_table_delimiter, 1)[0]
        feature = feature.split(config.gene_table_category_delimiter, 1)[0]
        feature = feature.split(config.name_mapping_join, 1)[0].rstrip("\r\n")
        if not self.warned and feature in config.legacy_pathway_feature_names:
            self.warned = True
            logger.warning(
                "DEPRECATION: This file (%s) seems to have been created with an earlier "
                "version of HUMAnN. %s and %s fractions in pathway-level data are being "
                "deprecated in favor of CPM-scale abundances emitted by the main workflow. "
                "This warning does not interrupt execution, and legacy input values remain "
                "supported. Users are advised to move away from previous pathway "
                "calculations in favor of CPM-based analyses.",
                self.path, *config.legacy_pathway_feature_names)


def warn_legacy_profile(lines, path):
    """Inspect a profile stream while preserving every input line."""
    warning = LegacyProfileWarning(path)
    for line in lines:
        warning.check(line)
        yield line
