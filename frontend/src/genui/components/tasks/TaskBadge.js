import { Badge } from 'reactstrap';
import React from 'react';

function TaskBadge(props) {
  if (props.tasks.length === 0) return null;

  return (
    <Badge {...props}>{props.tasks.length} {props.children}</Badge>
  )
}

export default TaskBadge;