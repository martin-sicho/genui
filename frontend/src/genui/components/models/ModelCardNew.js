import React from 'react';
import { CardHeader } from 'reactstrap';
import ModelFormRenderer from './ModelFormRenderer';
import ModelFormCardBody from './ModelFormCardBody';

class ModelCardNew extends React.Component {

  newModelFromFormData = (data) => {
    fetch(
      this.props.listURL
      , {
        method: 'POST'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    ).then((data) => this.props.handleResponseErrors(data, "Request failed. Data wrong or incomplete?")).then(
      data => {
        this.props.handleAddModel(this.props.modelClass, data)
      }
    ).catch(
      error => console.log(error)
    );
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>Create New {this.props.chosenAlgorithm.name} Model</CardHeader>
        <ModelFormRenderer
          {...this.props}
          formComponent={ModelFormCardBody}
          handleCreate={this.newModelFromFormData}
          project={this.props.currentProject}
        />
      </React.Fragment>
    )
  }
}

export default ModelCardNew;