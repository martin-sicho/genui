# Generated by Django 2.2.8 on 2020-01-29 16:04

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('generators', '0001_initial'),
    ]

    operations = [
        migrations.AddField(
            model_name='drugexnet',
            name='parent',
            field=models.ForeignKey(null=True, on_delete=django.db.models.deletion.CASCADE, to='generators.DrugExNet'),
        ),
        migrations.AlterField(
            model_name='drugeexcorpus',
            name='network',
            field=models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='corpus_set', to='generators.DrugExNet'),
        ),
    ]
