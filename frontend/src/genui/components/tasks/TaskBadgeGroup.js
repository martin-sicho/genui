import React from 'react';
import TaskBadge from './TaskBadge';

function TaskBadgeGroup(props) {
  const tasks = props.tasks;
  return (
    <React.Fragment>
      <TaskBadge href="#" color="primary" tasks={tasks.running}>Running</TaskBadge> <TaskBadge href="#" color="success" tasks={tasks.completed}>Completed</TaskBadge> <TaskBadge href="#" color="danger" tasks={tasks.errors}>Failed</TaskBadge>
    </React.Fragment>
  )
}

export default TaskBadgeGroup;