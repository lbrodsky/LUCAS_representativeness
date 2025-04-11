import os
import logging

def latest_version(dirs):
    v_max = -1
    for cntr in dirs:
        if not cntr.is_dir():
            continue
        versions = cntr.glob('v*')
        v_max = 0
        for v in versions:
            v_num = int(v.name[1:])
            v_max = max(v_num, v_max)

    return v_max

country_codes = {
    "austria": "AT",
    "belgium": "BE",
    "bulgaria": "BG",
    "croatia": "HR",
    "cyprus": "CY",
    "czech-republic": "CZ",
    "denmark": "DK",
    "estonia": "EE",
    "finland": "FI",
    "france": "FR",
    "germany": "DE",
    "great-britain": "UK",
    "greece": "EL",
    "hungary": "HU",
    "ireland": "IE",
    "italy": "IT",
    "latvia": "LV",
    "lithuania": "LT",
    "luxembourg": "LU",
    "malta": "MT",
    "netherlands": "NL",
    "poland": "PL",
    "portugal": "PT",
    "romania": "RO",
    "slovakia": "SK",
    "slovenia": "SI",
    "spain": "ES",
    "sweden": "SE"
}
