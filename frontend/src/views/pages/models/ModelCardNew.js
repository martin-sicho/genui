import React from "react"
import { CardHeader } from 'reactstrap';
import ModelCreateForm from './CreateForm';

class ModelCardNew extends React.Component {

  newModelFromFormData = (data) => {
    // TODO: handle this (do not forget to pass the new model up to the object handler via "handleAddModel")
    console.log(data)
  };

  render() {
    let molsets = [];
    Object.keys(this.props.compoundSets).forEach(
      (key) => molsets = molsets.concat(this.props.compoundSets[key])
    );

    if (molsets.length === 0) {
      // TODO: do something smarter like display an alert or redirect to the compounds UI
      return <p>There are no existing compound sets. Create one first.</p>
    }

    return (
      <React.Fragment>
        <CardHeader>Create New {this.props.chosenAlgorithm.name}</CardHeader>
        <ModelCreateForm
          handleCreate={this.newModelFromFormData}
          chosenAlgorithm={this.props.chosenAlgorithm}
          molsets={molsets}
          descriptors={this.props.descriptors}
        />
      </React.Fragment>
    )
  }
}

export default ModelCardNew;