from rest_framework import serializers
from django.contrib.auth.models import User
from .models import Hormone, Literature, AnalysisHistory, AnalysisTask


class HormoneSerializer(serializers.ModelSerializer):
    class Meta:
        model = Hormone
        fields = '__all__'


class LiteratureSerializer(serializers.ModelSerializer):
    class Meta:
        model = Literature
        fields = '__all__'


class AnalysisHistorySerializer(serializers.ModelSerializer):
    class Meta:
        model = AnalysisHistory
        fields = '__all__'


class UserSerializer(serializers.ModelSerializer):
    password = serializers.CharField(write_only=True)
    
    class Meta:
        model = User
        fields = ['id', 'username', 'email', 'password']
    
    def create(self, validated_data):
        user = User.objects.create_user(
            username=validated_data['username'],
            email=validated_data.get('email', ''),
            password=validated_data['password']
        )
        return user    