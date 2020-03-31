import React from 'react';
import { Progress } from 'reactstrap';

class TaskProgressBar extends React.Component {
  constructor(props) {
    super(props);

    this.progressURL =this.props.progressURL;

    this.state = {
      progressData : []
    }
  }

  componentDidMount() {
    this.updateProgress();
    this.interval = setInterval(this.updateProgress, 2000);
  }

  componentWillUnmount() {
    clearInterval(this.interval);
  }

  updateProgress = () => {
    const progressData = [];
    this.props.tasks.forEach(task => {
      const url = new URL(task.task_id + '/', this.progressURL);
      fetch(url, {credentials: "include",})
        .then(response => response.json())
        .then(data => {
          data.task = task;
          progressData.push(data);

          // FIXME: this set state should not happen if the component is unmounted -> the fetch needs to be cancelled properly
          this.setState(state => {
            return {
              progressData : progressData
            };
          });
        })
      ;
    });
  };

  render() {
    const tasks = this.props.tasks;
    if (tasks.length === 0) {
      return null;
    }

    const progress = this.state.progressData;
    progress.sort((a, b) => (a.task.task_id > b.task.task_id) ? 1 : -1);
    return (
      <React.Fragment>
        {
          progress.map(data => (
            <React.Fragment key={data.task.task_id}>
              <div className="text-center">{data.task.task_name} ({data.progress.percent}%)</div>
              <Progress value={data.progress.percent} />
            </React.Fragment>
          ))
        }
      </React.Fragment>
    )
  }
}

export default TaskProgressBar;