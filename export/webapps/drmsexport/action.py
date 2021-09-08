#!/usr/bin/env python3

from enum import Enum

__all__ = [ 'ErrorCode', 'ActionError', 'ActionApiError', 'Action' ]

class ErrorCode(Enum):
    UNKNOWN = (1, 'unknown error')
    API_ERROR = (2, 'API error')

    def __new__(cls, value, name):
        member = object.__new__(cls)
        member._value_ = value
        member.fullname = name
        return member

    def __int__(self):
        return self.value

class ActionError(Exception):
    _error_code = ErrorCode.UNKNOWN

    def __init__(self, *, msg='generic error'):
        self._msg = msg

    def __str__(self):
        return self._msg

class ActionApiError(ActionError):
    _error_code = ErrorCode.API_ERROR

    def __init__(self, *, msg=None):
        super().__init__(msg=msg)

class ActionType(type):
    _required_attributes = [ 'actions' ]
    _methods = set()
    _classes = {}

    def __new__(cls, name, bases, attributes):
        for required_att in cls._required_attributes:
            if required_att not in attributes:
                raise ActionApiError(msg=f'missing required attribute {required_att} in ActionType class {str(cls)}')

        for action in attributes['actions']:
            if action in cls._methods:
                raise ActionApiError(msg=f'duplicate action {action} defined in ActionType class {str(cls)}')
            elif action not in attributes:
                raise ActionApiError(msg=f'action {action} not defined in ActionType class {str(cls)}')
            else:
                cls._methods.add(action)

        cls = super(ActionType, cls).__new__(cls, name, bases, attributes)

        return cls

    def __init__(cls, name, bases, attributes):
        for action in attributes['actions']:
            cls._classes[action] = cls

class Action(metaclass=ActionType):
    # abstract class handles no actions - tuple is immutable
    actions = ()

    def __call__(self):
        return self._method()

    @classmethod
    def action(cls, *, action_type, args):
        action = cls._classes[action_type](method=action_type, **args)
        return action

    @classmethod
    def perform_action(cls, *, action_type, args):
        action = cls._classes[action_type](method=action_type, **args)
        return action.__call__()
