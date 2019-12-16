import React from 'react';
import {Redirect} from "react-router-dom";

class ProjectOverview extends React.Component {

    constructor(props) {
        super(props);
        this.state = {
            notFound : false
            , project : props.currentProject
        };
    }

    componentDidMount() {
        const { name } = this.props.match.params;
        const url = new URL(name + '/', this.props.apiUrls.projectList);
        fetch(url)
        .then(response => response.json())
        .then(this.processProjectData);
    }

    processProjectData = (data) => {
        const project = Object.assign({url : this.props.location.pathname}, data);
        if (project.id) {
          document.title = project.name;
          this.props.onProjectOpen(project);
          this.setState({project : project});
        } else {
          this.setState({notFound : true});
        }
    };

    render() {
        if (this.state.notFound) {
            return <Redirect to='/404' />
        }

        if (this.state.project) {
            return(
                this.notFound ? <Redirect to='/404' /> :
                <div>
                  <p>{`${this.state.project.name} overview...`}</p>
                </div>
            )
        } else {
            return <div>Loading...</div>
        }
    };
}

export default ProjectOverview;