#!/usr/bin/env python3

from enum import Enum

class StatusCode(Enum):
    UNKNOWN = 1, 'unknown'

    def __new__(cls, value, name):
        member = object.__new__(cls)
        member._value_ = value
        member.fullname = name
        return member

    def __int__(self):
        return self.value
