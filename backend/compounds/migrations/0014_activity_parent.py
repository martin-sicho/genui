# Generated by Django 2.2.8 on 2020-04-10 13:44

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0013_auto_20200409_1728'),
    ]

    operations = [
        migrations.AddField(
            model_name='activity',
            name='parent',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='children', to='compounds.Activity'),
        ),
    ]