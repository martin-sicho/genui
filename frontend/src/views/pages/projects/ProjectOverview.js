import React from 'react';

class ProjectOverview extends React.Component {

    render() {
        if (this.props.currentProject) {
            document.title = this.props.currentProject.name;
            return(
                <div>
                  <h1>{this.props.currentProject.name}</h1>
                  <p>
                      This page is not finished, but it will contain the summary of objects and tasks within the project in the future.
                  </p>
                </div>
            )
        } else {
            return <div>Loading...</div>
        }
    };
}

export default ProjectOverview;