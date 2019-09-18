from os import environ as _env

_off = "scalar_off"

if _off in _env: del _env[_off]

#_env[_off] = "yes" # for testing

from units import *
