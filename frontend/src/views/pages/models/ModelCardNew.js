import React from "react"
import { CardHeader } from 'reactstrap';
import ModelCreateForm from './CreateForm';

class ModelCardNew extends React.Component {

  newModelFromFormData = (data) => {
    // TODO: handle this (do not forget to pass the new model up to the object handler via "handleAddModel")
    console.log(data)
  };

  render() {
    return (
      <React.Fragment>
        <CardHeader>Create New {this.props.newModel.name}</CardHeader>
        <ModelCreateForm
          handleCreate={this.newModelFromFormData}
          newModel={this.props.newModel}
          molsets={this.props.compoundSets}
        />
      </React.Fragment>
    )
  }
}

export default ModelCardNew;