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
    percent = serializers.IntegerField()

class TaskProgressSerializer(serializers.Serializer):
    complete = serializers.BooleanField()
    success = serializers.BooleanField()
    progress = TaskProgressInfoSerializer()
    result = serializers.DictField()

