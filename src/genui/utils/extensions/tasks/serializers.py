"""
serializers

Created by: Martin Sicho
On: 4/30/20, 5:33 PM
"""
from rest_framework import serializers


class TaskSerializer(serializers.Serializer):
    task_id = serializers.CharField(allow_blank=False)
    status = serializers.CharField(allow_blank=False)
    result = serializers.JSONField(allow_null=True)
    traceback = serializers.CharField(allow_null=True, max_length=100000)


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