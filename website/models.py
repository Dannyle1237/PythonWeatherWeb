from . import db 
from flask_login import UserMixin
from sqlalchemy.sql import func

class City(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(10000))
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))

class User(db.Model, UserMixin):
    id = db.Column(db.Integer, primary_key=True)
    email = db.Column(db.String(150), unique=True)
    password = db.Column(db.String(150))
    first_name = db.Column(db.String(150))
    cities = db.relationship('City')
    weather = db.relationship('CityWeather')

class CityWeather(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    name = db.Column(db.String(10000))
    temp_k = db.Column(db.Integer)
    temp_f = db.Column(db.Numeric(precision=10, scale=2))
    temp_c = db.Column(db.Numeric(precision=10, scale=2))
    temp_min = db.Column(db.Numeric(precision=10, scale=2))
    temp_max = db.Column(db.Numeric(precision=10, scale=2))
    description = db.Column(db.String(150))
    wind_speed = db.Column(db.Numeric(precision=10, scale=2))
    wind_dir = db.Column(db.Integer)
    user_id = db.Column(db.Integer, db.ForeignKey('user.id'))
