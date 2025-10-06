"""Load configuration file for pypima."""

import os
from configparser import ConfigParser
from pathlib import Path


class QuoteStrippingConfigParser(ConfigParser):
    """ConfigParser subclass that strips quotes from option values."""

    def get(self, section: str, option: str, *args, **kwargs) -> str:
        """Get an option value for a given section stripping quotes."""
        value = super().get(section, option, *args, **kwargs)
        return value.strip('"')

    def getpath(self, section: str, option: str, *args, **kwargs) -> Path:
        """Get an option value for a given section as a Path object."""
        value = super().get(section, option, *args, **kwargs)
        return Path(value.strip('"')).expanduser()


def load_config(
    file_path: str | os.PathLike | None = None,
) -> QuoteStrippingConfigParser:
    """Load configuration file."""
    config = QuoteStrippingConfigParser()

    if file_path is None:
        config_dir = Path(os.getenv("XDG_CONFIG_HOME", "~/.config")).expanduser()
        file_path = config_dir / "pypima.cfg"

    with open(file_path) as file:
        config.read_file(file)

    return config
