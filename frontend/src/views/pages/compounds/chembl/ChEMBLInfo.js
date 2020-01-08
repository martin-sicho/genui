import React from 'react';
import { Badge, Col, Progress, Row } from 'reactstrap';

class MolsStats extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      molCount : 0,
      lastUpdate : null
    }
  }

  componentDidMount() {
    this.fetchUpdates();
  }

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.props.molsetIsUpdating) {
      this.fetchUpdates();
    }
  }

  fetchUpdates = () => {
    const time_now = new Date().getTime();
    if (!this.state.lastUpdate || (time_now - this.state.lastUpdate > 2000)) {
      // TODO: add proper request error handling
      fetch(this.props.moleculesURL)
        .then(response => response.json())
        .then(this.updateMolStats)
    }
  };

  updateMolStats = (data) => {
    this.setState({
      molCount : data.count
      , lastUpdate : new Date().getTime()
    });
  };

  render() {
    return (
      <React.Fragment>
        <h4>Compounds</h4>
        <p>Unique in Total: {this.state.molCount}</p>
        <h4>Associated Targets</h4>
        <ul>
          {
            this.props.molset.targets.map(
              target => <li key={target.targetID}>{target.targetID}</li>
            )
          }
        </ul>
      </React.Fragment>
    )
  }
}

function TaskHeadingBadge(props) {
  if (props.tasks.length === 0) return null;

  return (
    <Badge {...props}>{props.tasks.length} {props.children}</Badge>
  )
}

class TasksProgressOverview extends React.Component {
  constructor(props) {
    super(props);

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
    const tasks = this.props.tasks;
    const progressData = [];
    tasks.forEach(task => {
      const url = new URL(task.task_id + '/', this.props.progressURL);
      fetch(url)
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

class MolSetTasksStatus extends React.Component {
  intervalID;

  constructor(props) {
    super(props);

    this.state = {
      tasks : null
    }
  }

  componentWillUnmount() {
    clearTimeout(this.intervalID);
  }

  componentDidMount() {
    this.updateTasks();
  }

  updateTasks = () => {
    fetch(this.props.tasksURL)
      .then(response => this.props.handleResponseErrors(response, 'Failed to fetch task info from backend.'))
      .then(data => {
        const tasks = this.groupTasks(data);
        this.setState({tasks : tasks});
        this.intervalID = setTimeout(this.updateTasks, 5000);
        this.props.processTasks(tasks);
      }).catch(
      (error) => console.log(error)
    )
  };

  groupTasks = (data) => {
    const tasks = data;
    const completed = [];
    const running = [];
    const errors = [];
    Object.keys(tasks).forEach(task_name => {
      tasks[task_name].forEach(task => {
        task.task_name = task_name;
        if (task.status === 'SUCCESS') {
          completed.push(task)
        } else if (['STARTED', 'RECEIVED', 'PENDING', 'RETRY', 'PROGRESS'].includes(task.status)) {
          running.push(task)
        } else if (['FAILURE', 'REVOKED'].includes(task.status)) {
          errors.push(task)
        }
      });
    });

    return {
      completed : completed,
      running : running,
      errors : errors
    }
  };

  render() {
    const tasks = this.state.tasks;
    if (!tasks) {
      return null
    }

    return (
      <React.Fragment>
        <h4>
          Tasks <TaskHeadingBadge href="#" color="primary" tasks={tasks.running}>Running</TaskHeadingBadge> <TaskHeadingBadge href="#" color="success" tasks={tasks.completed}>Completed</TaskHeadingBadge> <TaskHeadingBadge href="#" color="danger" tasks={tasks.errors}>Failed</TaskHeadingBadge>
        </h4>
        {
          tasks.running.length > 0 ?
          <TasksProgressOverview {...this.props} tasks={tasks.running}/>
          : null
        }
      </React.Fragment>
    )
  }
}

class ChEMBLInfo extends React.Component {

  constructor(props) {
    super(props);

    this.molset = this.props.molset;
  }

  render() {
    return (
      <Row>
        <Col sm="12">
          <h4>Description</h4>
          <p>{this.molset.description}</p>
          <MolsStats molset={this.props.molset} moleculesURL={this.props.moleculesURL} molsetIsUpdating={this.props.molsetIsUpdating} />
          <MolSetTasksStatus
            progressURL={this.props.apiUrls.celeryProgress}
            tasksURL={this.props.tasksURL}
            molset={this.props.molset}
            processTasks={this.props.processTasks}
            handleResponseErrors={this.props.handleResponseErrors}
            molsetIsUpdating={this.props.molsetIsUpdating}
          />
        </Col>
      </Row>
    );
  }
}

export default ChEMBLInfo;