from aiida.orm import Data
import datetime


class DateTimeData(Data):
    """AiiDA node to store a datetime.datetime object."""

    def __init__(self, value: datetime.datetime, **kwargs):
        if not isinstance(value, datetime.datetime):
            raise TypeError(f"Expected datetime.datetime, got {type(value)}")
        super().__init__(**kwargs)
        # Store as ISO string for portability
        self.base.attributes.set("datetime", value.isoformat())

    @property
    def value(self) -> datetime.datetime:
        """Return the stored datetime as a datetime object."""
        return datetime.datetime.fromisoformat(self.base.attributes.get("datetime"))

    def __str__(self):
        return str(self.value)
