"""
serializers

Created by: Martin Sicho
On: 12/22/19, 6:26 PM
"""
from rest_framework import serializers


class TaskSerializer(serializers.Serializer):
    task_id = serializers.CharField(allow_blank=False)
    status = serializers.CharField(allow_blank=False)

class TasksSerializerFactory:
    class AutoSchemaMixIn:
        def get_operation(self, path, method):
            ret = super().get_operation(path, method)
            ret['responses']['200']['content']['application/json']['schema'] = {
                'type' : 'object',
                'properties' : {
                        'taskName' : {
                            'type' : 'array',
                            'items' : {
                                'properties' : {
                                    'task_id' : {
                                        'type' : 'string'
                                    },
                                    'status' : {
                                        'type' : 'string'
                                    },
                                    'required': ['task_id', 'status']
                                }
                            }
                        },
                        'required': ['taskName']
                    }
                }
            return ret

    @staticmethod
    def get(field_names):
        return type('TaskSetSerializer', (serializers.Serializer,), {x : serializers.ListField(child=TaskSerializer(required=False), allow_empty=True) for x in field_names})

class TaskProgressInfoSerializer(serializers.Serializer):
    current = serializers.IntegerField()
    total = serializers.IntegerField()
    percent = serializers.FloatField()

class TaskProgressSerializer(serializers.Serializer):
    success = serializers.BooleanField(allow_null=True)
    complete = serializers.BooleanField()
    progress = TaskProgressInfoSerializer()

class GenericModelSerializerMixIn:
    className = serializers.CharField(default="")
    extraArgs = serializers.DictField(required=False, allow_empty=True)

    def to_representation(self, instance):
        ret = super().to_representation(instance)
        model_class = instance.__class__
        ret['className'] = model_class.__name__
        base_fields = set(self.Meta.model._meta.get_fields())
        derived_fields = set(instance.__class__._meta.get_fields())
        extra_fields = []
        for x in derived_fields - base_fields:
            if not x.name.endswith("_ptr"):
                extra_fields.append(x.name)

        if extra_fields:
            serializer_class = type(
                'GenericModelSerializer'
                , (serializers.Serializer,),
                {x : serializers.ModelField(model_field=model_class._meta.get_field(x)) for x in extra_fields})
            extra_data = serializer_class(instance).data
            ret['extraArgs'] = extra_data
        else:
            ret['extraArgs'] = {}
        return ret