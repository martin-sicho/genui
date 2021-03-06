# Generated by Django 2.2.8 on 2020-05-05 07:27

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('models', '0001_initial'),
    ]

    operations = [
        migrations.AlterField(
            model_name='model',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_models.model_set+', to='contenttypes.ContentType'),
        ),
        migrations.AlterField(
            model_name='modelparametervalue',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_models.modelparametervalue_set+', to='contenttypes.ContentType'),
        ),
        migrations.AlterField(
            model_name='modelperformance',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_models.modelperformance_set+', to='contenttypes.ContentType'),
        ),
        migrations.AlterField(
            model_name='trainingstrategy',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_models.trainingstrategy_set+', to='contenttypes.ContentType'),
        ),
        migrations.AlterField(
            model_name='validationstrategy',
            name='polymorphic_ctype',
            field=models.ForeignKey(editable=False, null=True, on_delete=django.db.models.deletion.CASCADE, related_name='polymorphic_models.validationstrategy_set+', to='contenttypes.ContentType'),
        ),
    ]
