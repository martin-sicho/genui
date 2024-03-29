# Generated by Django 2.2.20 on 2021-04-28 14:26

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0009_auto_20200728_1326'),
    ]

    operations = [
        migrations.CreateModel(
            name='MolSetExporter',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=16, unique=True)),
                ('classPath', models.CharField(max_length=512, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='MolSetExport',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=128)),
                ('description', models.TextField(blank=True, max_length=10000)),
                ('exporter', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, to='compounds.MolSetExporter')),
                ('molset', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='exports', to='compounds.MolSet')),
            ],
        ),
        migrations.AddField(
            model_name='molsetfile',
            name='export',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='files', to='compounds.MolSetExport'),
        ),
    ]
