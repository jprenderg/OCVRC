"""Shared query-management settings for optical/IR/UV/radio catalogs."""

from _query_manager import run_queries

def managed(survey, sources, query_one, output, radius, section, config, start, end):
    return run_queries(
        survey, sources, query_one, output, radius=radius,
        retries=config.getint(section, 'MAX_RETRIES', fallback=2),
        delay=config.getfloat(section, 'RETRY_DELAY_SECONDS', fallback=3),
        failure_limit=config.getint('QUERY_MANAGEMENT', 'CONSECUTIVE_FAILURE_LIMIT', fallback=5),
        save_every=config.getint('QUERY_MANAGEMENT', 'SAVE_EVERY', fallback=10),
        start=max(0, start), end=end,
        resume=config.getboolean('QUERY_MANAGEMENT', 'RESUME', fallback=True))
