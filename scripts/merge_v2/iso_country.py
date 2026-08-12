#!/usr/bin/env python3
"""ISO 3166-1 alpha-2 → (标准国名, 洲)。pycountry 未装，静态表（覆盖 MBM 68 码 + 常见）。
洲名对齐 canonical：Europe/North America/Asia/Africa/Oceania/South America。"""
import pandas as pd

# code -> (country, continent)
ISO2 = {
    "US":("United States","North America"), "CN":("China","Asia"), "DK":("Denmark","Europe"),
    "JP":("Japan","Asia"), "CA":("Canada","North America"), "GB":("United Kingdom","Europe"),
    "FI":("Finland","Europe"), "DE":("Germany","Europe"), "BD":("Bangladesh","Asia"),
    "IN":("India","Asia"), "NL":("Netherlands","Europe"), "SE":("Sweden","Europe"),
    "RU":("Russia","Europe"), "NZ":("New Zealand","Oceania"), "IT":("Italy","Europe"),
    "ES":("Spain","Europe"), "AU":("Australia","Oceania"), "NG":("Nigeria","Africa"),
    "BE":("Belgium","Europe"), "AT":("Austria","Europe"), "FR":("France","Europe"),
    "IL":("Israel","Asia"), "MW":("Malawi","Africa"), "TH":("Thailand","Asia"),
    "IE":("Ireland","Europe"), "CM":("Cameroon","Africa"), "NO":("Norway","Europe"),
    "TW":("Taiwan","Asia"), "HK":("Hong Kong","Asia"), "MX":("Mexico","North America"),
    "SG":("Singapore","Asia"), "MG":("Madagascar","Africa"), "TZ":("Tanzania","Africa"),
    "KR":("South Korea","Asia"), "CH":("Switzerland","Europe"),
    "CF":("Central African Republic","Africa"), "RO":("Romania","Europe"),
    "EE":("Estonia","Europe"), "MZ":("Mozambique","Africa"), "ZW":("Zimbabwe","Africa"),
    "UG":("Uganda","Africa"), "PL":("Poland","Europe"), "CO":("Colombia","South America"),
    "GR":("Greece","Europe"), "KE":("Kenya","Africa"), "BR":("Brazil","South America"),
    "EC":("Ecuador","South America"), "TR":("Turkey","Asia"), "GH":("Ghana","Africa"),
    "VE":("Venezuela","South America"), "PT":("Portugal","Europe"), "ML":("Mali","Africa"),
    "KZ":("Kazakhstan","Asia"), "GT":("Guatemala","North America"), "RS":("Serbia","Europe"),
    "SD":("Sudan","Africa"), "BW":("Botswana","Africa"), "AZ":("Azerbaijan","Asia"),
    "JO":("Jordan","Asia"), "LT":("Lithuania","Europe"), "HU":("Hungary","Europe"),
    "UM":("United States Minor Outlying Islands","Oceania"), "MY":("Malaysia","Asia"),
    "HR":("Croatia","Europe"), "PE":("Peru","South America"), "IS":("Iceland","Europe"),
    "AM":("Armenia","Asia"), "AE":("United Arab Emirates","Asia"),
    # 常见补充
    "PK":("Pakistan","Asia"), "ID":("Indonesia","Asia"), "VN":("Vietnam","Asia"),
    "PH":("Philippines","Asia"), "EG":("Egypt","Africa"), "ZA":("South Africa","Africa"),
    "AR":("Argentina","South America"), "CL":("Chile","South America"), "UA":("Ukraine","Europe"),
    "CZ":("Czechia","Europe"), "SK":("Slovakia","Europe"), "BG":("Bulgaria","Europe"),
    "SI":("Slovenia","Europe"), "LV":("Latvia","Europe"), "LU":("Luxembourg","Europe"),
    "SA":("Saudi Arabia","Asia"), "IR":("Iran","Asia"), "IQ":("Iraq","Asia"),
    "ET":("Ethiopia","Africa"), "GM":("Gambia","Africa"), "SN":("Senegal","Africa"),
    "MA":("Morocco","Africa"), "BF":("Burkina Faso","Africa"), "CI":("Côte d'Ivoire","Africa"),
    "FJ":("Fiji","Oceania"), "PG":("Papua New Guinea","Oceania"),
}

def country_of(iso):
    if pd.isna(iso): return pd.NA
    return ISO2.get(str(iso).strip().upper(), (pd.NA, pd.NA))[0]

def continent_of(iso):
    if pd.isna(iso): return pd.NA
    return ISO2.get(str(iso).strip().upper(), (pd.NA, pd.NA))[1]
