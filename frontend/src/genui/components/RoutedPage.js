import React from "react"
import {Redirect} from "react-router-dom";
import PageAlertContext from '../../vibe/components/PageAlert/PageAlertContext';
import '../styles.css'

/*
 * Component which serves the purpose of a "root route component".
 * 
 * Source: https://stackoverflow.com/a/54112771
 */
class RoutedPage extends React.Component {

  constructor(props) {
    super(props);
    this.state = {
      notFound : false
      , project : null
    };
  }

  sleep = (milliseconds) => {
    return new Promise(resolve => setTimeout(resolve, milliseconds))
  };

  handleResponseErrors = (
    response,
    message='Failed to fetch data from backend.',
    showAlert=false,
    callback=null
  ) => {
    if (!response.ok) {
      if (showAlert) {
        this.showAlert(message);
      }
      console.log(response);
      response.json()
        .then(data => {
          console.log(data);
          if (callback) {
            callback(data);
          }
        })
        .catch(e => console.log(e));
      throw new Error(message);
    } else {
      return response.json();
    }
  };

  showAlert = (message, severity='danger') => {
    this.context.setAlert(message, severity);
  };

  retryAction = (action, message='', severity='danger', interval=5000) => {
    if (message) {
      this.showAlert(message + ` (retrying in ${interval / 1000} seconds)`, severity);
    }
    console.log(message + ` Retrying in ${interval / 1000} seconds...`);
    this.sleep(interval)
      .then(action)
    ;
  };

  setPageTitle = (newTitle) => {
    this.props.handlePageTitleChange(newTitle)
  };
  
  /**
   * Here, we define a react lifecycle method that gets executed each time
   * our component is mounted to the DOM, which is exactly what we want in this case
   */
  componentDidMount() {
    this.fetchProject();
    this.setPageTitle(this.props.title);
  }

  componentWillUnmount() {
    this.props.onHeaderChange(null);
  }

  fetchProject = () => {
    const { project } = this.props.match.params;
    if (project) {
      const url = new URL(project + '/', this.props.apiUrls.projectList);
      fetch(url, {credentials: 'include'})
        .then((response) => this.handleResponseErrors(response, 'Failed to fetch project data from backend.'))
        .then(this.processProjectData)
        .catch(
          () => {
            this.retryAction(this.fetchProject, 'Failed to fetch project data from backend.')
          }
        );
    }
  };
  
   processProjectData = (data) => {
        const project = Object.assign({url : `/projects/${data.id}/`}, data);
        if (project.id) {
          this.props.onProjectOpen(project);
          this.setState({project : project});
        } else {
          this.setState({notFound : true});
        }
    };

  /**
   * Here, we use a component prop to render
   * a component, as specified in route configuration
   */
  render() {
    if (this.state.notFound) {
      return <Redirect to='/404' />
    }
    
    const PageComponent = this.props.component;
    return (
      <PageComponent
        {...this.props}
        currentProject={this.state.project}
        setPageTitle={this.setPageTitle}
        retryAction={this.retryAction}
        handleResponseErrors={this.handleResponseErrors}
      />
    )
  }
}
RoutedPage.contextType = PageAlertContext;

export default RoutedPage