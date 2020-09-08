# Generated by Django 2.2.8 on 2020-07-28 13:26

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0008_auto_20200723_1212'),
    ]

    operations = [
        migrations.AlterField(
            model_name='pictureformat',
            name='extension',
            field=models.CharField(choices=[('.png', 'PNG'), ('.svg', 'SVG')], default='.svg', max_length=128, unique=True),
        ),
    ]