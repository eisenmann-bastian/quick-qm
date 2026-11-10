"""Shared logging configuration for Quick-QM applications."""

import logging
import sys


LOGGER_NAME = "quick_qm"


def get_logger(name):
    """Return a logger under the shared Quick-QM namespace."""
    return logging.getLogger(f"{LOGGER_NAME}.{name}")


def configure_logging(level=logging.INFO, stream=None):
    """Configure console logging once, without changing unrelated loggers."""
    logger = logging.getLogger(LOGGER_NAME)
    if not any(getattr(handler, "_quick_qm_handler", False) for handler in logger.handlers):
        handler = logging.StreamHandler(stream if stream is not None else sys.stderr)
        handler.setFormatter(logging.Formatter("%(levelname)s %(name)s: %(message)s"))
        handler._quick_qm_handler = True
        logger.addHandler(handler)
    logger.setLevel(level)
    logger.propagate = False
