import datetime as dt
import requests
import os

BASE_URL="http://api.openweathermap.org/data/2.5/weather?"
API_KEY=open(os.getcwd() + '/website/api/apikey.txt', 'r').read()

def requestItemMain(url, item):
    return requests.get(url).json()['main'][item]


def check_if_city_exists(CITY):
    url = BASE_URL + "appid=" + API_KEY + "&q=" + CITY
    return requests.get(url)
    
def kelvin_to_celsius_fahrenheit(kelvin):
    celsius = kelvin-273.15
    fahrenheit = celsius*(9/5) + 32
    return celsius, fahrenheit, kelvin

#Return temps in order C, F, K, Max, Min
def getTemps(CITY):
    url = BASE_URL + "appid=" + API_KEY + "&q=" + CITY
    temp_kelvin = requestItemMain(url, "temp")
    temps = kelvin_to_celsius_fahrenheit(temp_kelvin)
    temp_max = requestItemMain(url, 'temp_max')
    temp_min = requestItemMain(url, 'temp_min')
    return temps, temp_max, temp_min

#Return description of sky
def getDescription(CITY):
    url = BASE_URL + "appid=" + API_KEY + "&q=" + CITY
    return requests.get(url).json()['weather'][0]['description']

#Wind speed and degree
def getWind(CITY):
    url = BASE_URL + "appid=" + API_KEY + "&q=" + CITY
    wind = requests.get(url).json()['wind']
    speed = wind['speed']
    dir = wind['deg']
    return speed, dir

