import React from 'react';

class ProjectOverview extends React.Component {

    render() {
        if (this.props.currentProject) {
            document.title = this.props.currentProject.name;
            return(
                <div>
                  <p>{`${this.props.currentProject.name} overview...`}</p>
                </div>
            )
        } else {
            return <div>Loading...</div>
        }
    };
}

export default ProjectOverview;