"""Logging utilities for finder pipeline scripts."""

from __future__ import annotations

import logging
import os
import traceback
from datetime import datetime
from enum import IntEnum
from functools import wraps
from typing import Any, Callable, Optional

from rich.progress import (
    BarColumn,
    Progress,
    TextColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)


class Logger:
    """Configure console and file logging with optional Rich progress bars.

    Attributes:
        LOG_DIR: Default directory for log files, relative to the working directory.
    """

    LOG_DIR = "logs"

    class Level(IntEnum):
        """Logging severity levels compatible with the standard logging module."""

        NOTSET = logging.NOTSET
        DEBUG = logging.DEBUG
        INFO = logging.INFO
        WARNING = logging.WARNING
        ERROR = logging.ERROR
        CRITICAL = logging.CRITICAL

    def __init__(
        self,
        name: Optional[str] = None,
        level: Level = Level.INFO,
        log_to_console: bool = True,
        log_to_file: bool = True,
        draw_progress: bool = False,
    ) -> None:
        """Initialize a named logger with optional console, file, and progress output.

        Args:
            name: Logger name passed to ``logging.getLogger``.
            level: Minimum message severity to emit.
            log_to_console: When True, attach a stream handler.
            log_to_file: When True, attach a timestamped file handler.
            draw_progress: When True, start a Rich progress display.
        """
        self.logger = logging.getLogger(name)
        self.logger.setLevel(level)

        os.makedirs(self.LOG_DIR, exist_ok=True)
        log_filename = os.path.join(
            self.LOG_DIR,
            f"log_{datetime.now().strftime('%Y-%m-%d_%H-%M-%S')}.log",
        )

        formatter = logging.Formatter(
            "[%(asctime)s - %(levelname)s]: %(message)s",
            "%Y-%m-%d %H:%M:%S",
        )

        if not self.logger.handlers:
            if log_to_console:
                console_handler = logging.StreamHandler()
                console_handler.setFormatter(formatter)
                self.logger.addHandler(console_handler)

            if log_to_file:
                file_handler = logging.FileHandler(log_filename, mode="a", encoding="utf-8")
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)

        if draw_progress:
            self.progress = Progress(
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TimeElapsedColumn(),
                TimeRemainingColumn(),
            )
            self.progress.start()
        else:
            self.progress = None

    def get_logger(self) -> logging.Logger:
        """Return the configured ``logging.Logger`` instance.

        Returns:
            The underlying standard-library logger.
        """
        return self.logger

    def trace_call(self, func: Callable[..., Any]) -> Callable[..., Any]:
        """Return a wrapper that logs function entry, exit, and exceptions.

        Args:
            func: Callable to wrap.

        Returns:
            Wrapped callable with identical behavior and logging.
        """
        @wraps(func)
        def wrapper(*args: Any, **kwargs: Any) -> Any:
            func_name = func.__qualname__
            self.logger.info(f"[START] {func_name}")
            try:
                result = func(*args, **kwargs)
                self.logger.info(f"[END] {func_name}")
                return result
            except Exception as exc:
                self.logger.error(f"[ERROR in {func_name}]: {exc}")
                self.logger.debug(traceback.format_exc())
                raise

        return wrapper

    def log_error(self, message: str, exception: Optional[Exception] = None) -> None:
        """Log an error message and optional exception traceback.

        Args:
            message: Human-readable error description.
            exception: Optional exception to include at debug level.
        """
        if exception:
            self.logger.error(f"{message} | Exception: {exception}")
            self.logger.debug(traceback.format_exc())
        else:
            self.logger.error(message)

    def progress_task(self, description: str, total: int = 100) -> Optional[Any]:
        """Create a Rich progress task when progress display is enabled.

        Args:
            description: Short label shown in the progress bar.
            total: Expected number of progress steps.

        Returns:
            Rich task identifier, or None when progress display is disabled.
        """
        if self.progress:
            return self.progress.add_task(description, total=total)
        return None

    def advance_progress(self, task: Optional[Any], advance: int = 1) -> None:
        """Advance a Rich progress task by the given number of steps.

        Args:
            task: Task identifier returned by ``progress_task``.
            advance: Number of steps to add.
        """
        if self.progress and task is not None:
            self.progress.update(task, advance=advance)

    def finish_progress(self, task: Optional[Any]) -> None:
        """Mark a progress task complete and stop the progress display.

        Args:
            task: Task identifier returned by ``progress_task``.
        """
        if self.progress and task is not None:
            self.progress.update(task, completed=100)
            self.progress.stop()

    def resources_info(self, message: str) -> None:
        """Log a resource-related informational message.

        Provided for compatibility with Snakemake logging hooks.

        Args:
            message: Text to log at INFO level.
        """
        self.logger.info(message)

    def run_info(self, message: str) -> None:
        """Log a run-related informational message.

        Provided for compatibility with Snakemake logging hooks.

        Args:
            message: Text to log at INFO level.
        """
        self.logger.info(message)
