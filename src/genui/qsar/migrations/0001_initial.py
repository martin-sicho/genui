# Generated by Django 2.2.8 on 2020-04-27 14:28

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    initial = True

    dependencies = [
        ('compounds', '0001_initial'),
        ('modelling', '0001_initial'),
    ]

    operations = [
        migrations.CreateModel(
            name='DescriptorGroup',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('name', models.CharField(max_length=128, unique=True)),
            ],
        ),
        migrations.CreateModel(
            name='ModelActivity',
            fields=[
                ('activity_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='compounds.Activity')),
            ],
            options={
                'abstract': False,
                'base_manager_name': 'objects',
            },
            bases=('compounds.activity',),
        ),
        migrations.CreateModel(
            name='QSARTrainingStrategy',
            fields=[
                ('trainingstrategy_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='modelling.TrainingStrategy')),
                ('activityThreshold', models.FloatField(null=True)),
                ('activitySet', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.ActivitySet')),
                ('activityType', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.ActivityTypes')),
                ('descriptors', models.ManyToManyField(to='qsar.DescriptorGroup')),
            ],
            options={
                'abstract': False,
                'base_manager_name': 'objects',
            },
            bases=('modelling.trainingstrategy',),
        ),
        migrations.CreateModel(
            name='QSARModel',
            fields=[
                ('model_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='modelling.Model')),
                ('molset', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, related_name='models', to='compounds.MolSet')),
                ('predictionsType', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.ActivityTypes')),
                ('predictionsUnits', models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='compounds.ActivityUnits')),
            ],
            options={
                'abstract': False,
            },
            bases=('modelling.model',),
        ),
        migrations.CreateModel(
            name='ModelActivitySet',
            fields=[
                ('activityset_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='compounds.ActivitySet')),
                ('model', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='predictions', to='qsar.QSARModel')),
            ],
            options={
                'abstract': False,
            },
            bases=('compounds.activityset',),
        ),
    ]