#!/usr/bin/env python3

from statuscode import StatusCode
from response import ErrorResponse

__all__ = [ 'Error' ]

class Error(Exception):
    _error_code = StatusCode(StatusCode.UNKNOWN)

    def __init__(self, *, msg=None):
        self._msg = msg
        if not hasattr(self, '_header'):
            self._header = f'[ {self.__class__.__name__} ]'
        self._response = None

    def __str__(self):
        return self._msg

    def _generate_response(self):
        if self._msg is None:
            msg = f'{self._header}'
        else:
            msg=f'{self._header} {self._msg}'

        self._response = ErrorResponse(error_code=self._error_code, msg=msg)

    @property
    def response(self):
        if self._response is None:
            self._generate_response()

        return self._response
