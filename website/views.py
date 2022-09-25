from flask import Blueprint, render_template, flash, request, jsonify, redirect, url_for
from flask_login import login_required, current_user
from .models import City, CityWeather
from . import db
import json 
from .api import weatherAPI as wAPI

views = Blueprint('views', __name__)

@views.route('/', methods=['GET', 'POST'])
@login_required
def home():
    if request.method == 'POST':
        city = request.form.get('city')
        cityExists = City.query.filter_by(name=city).first()
        if cityExists:
            flash('You already have ' +city + ' listed', category = "error")
        elif not wAPI.check_if_city_exists(city):
            flash('City does not exist', category="eror")
        elif len(city) <= 1:
            flash('Please type in city name', category = "error")
        else:
            new_city = City(name=city, user_id=current_user.id)
            db.session.add(new_city)
            db.session.commit()
            flash('City added!', category='success')
    return render_template("home.html", user=current_user)

@views.route('/delete-note', methods=['POST'])
def delete_note():
    city = json.loads(request.data)
    cityId = city['cityId']
    city = City.query.get(cityId) 
    if city:
        if city.user_id == current_user.id:
            db.session.delete(city)
            db.session.commit()
    return jsonify({})

@views.route('/weather', methods=['GET','POST'])
def weather():
    if request.method == 'POST':
        weather = CityWeather.query.filter_by(user_id = current_user.id).first()
        if weather:
            print("Weather deleted")
            db.session.delete(weather)
            db.session.commit()
        city = json.loads(request.data)
        cityId = city['cityId']
        city_name = City.query.get(cityId).name
        temps = wAPI.getTemps(city_name)
        description = wAPI.getDescription(city_name)
        wind = wAPI.getWind(city_name)
        new_weather = CityWeather(name=city_name, temp_c=temps[0][0], temp_f=temps[0][1], temp_k=int(temps[0][2]), temp_max=temps[1], temp_min=temps[2], description=description, wind_speed=wind[0], wind_dir=int(wind[1]), user_id = current_user.id)
        db.session.add(new_weather)
        db.session.commit()

    return render_template("weather.html", user=current_user)