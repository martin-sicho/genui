# Generated by Django 2.2.8 on 2020-05-12 08:17

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0002_auto_20200427_2039'),
    ]

    operations = [
        migrations.AlterField(
            model_name='activityset',
            name='molecules',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='activities', to='compounds.MolSet'),
        ),
    ]
