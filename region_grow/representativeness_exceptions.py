#!/usr/bin/python3

from abc import ABC
import logging

class BaseError(Exception, ABC):
    """Error to be raised in the case of wrong parameter of a script config.
    """
    def __init__(self, message, *args, **kwargs):
        logging.error(f'{self.__class__.__name__}: {message}')

        super().__init__(message, *args, **kwargs)

class ConfigError(BaseError):
    """Error to be raised in the case of wrong parameter of a script config.
    """

    pass

class IllegalArgumentError(BaseError):
    """Error to be raised in the case of wrong argument in a function call.
    """

    pass
