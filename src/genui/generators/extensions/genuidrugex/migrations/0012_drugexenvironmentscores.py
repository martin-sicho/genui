# Generated by Django 4.1 on 2022-09-01 12:21

from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ('compounds', '0011_alter_activity_polymorphic_ctype_and_more'),
        ('genuidrugex', '0011_alter_drugexenvironment_rewardscheme'),
    ]

    operations = [
        migrations.CreateModel(
            name='DrugExEnvironmentScores',
            fields=[
                ('activityset_ptr', models.OneToOneField(auto_created=True, on_delete=django.db.models.deletion.CASCADE, parent_link=True, primary_key=True, serialize=False, to='compounds.activityset')),
                ('environment', models.ForeignKey(on_delete=django.db.models.deletion.CASCADE, related_name='scores', to='genuidrugex.drugexenvironment')),
            ],
            options={
                'abstract': False,
            },
            bases=('compounds.activityset',),
        ),
    ]