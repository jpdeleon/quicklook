"""Custom exceptions for the quicklook pipeline."""


class QuickLookError(Exception):
    """Base exception for all quicklook pipeline errors."""


class NoDataError(QuickLookError):
    """Raised when no data is found for a target (no lightcurves, TPFs, etc.)."""


class InvalidInputError(QuickLookError):
    """Raised when user-supplied input is invalid (bad ephemeris, missing params, etc.)."""


class PipelineError(QuickLookError):
    """Raised when the pipeline encounters an unrecoverable processing error."""
