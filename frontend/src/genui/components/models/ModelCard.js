import React from 'react';
import { Button, CardBody, CardFooter, CardHeader } from 'reactstrap';
import { ProjectItemSubTitle, TabWidget } from '../../index';

class ModelCard extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      isDeleting : false
    }
  }

  render() {
    const model = this.props.model;

    return (
      <React.Fragment>
        <CardHeader>{model.name}</CardHeader>

        <CardBody className="scrollable">
          <ProjectItemSubTitle item={model}/>
          <TabWidget {...this.props} tabs={this.props.tabs}/>
        </CardBody>

        <CardFooter>
          <Button color="danger" disabled={this.state.isDeleting} onClick={() => {
            this.setState({ isDeleting: true });
            this.props.onModelDelete(this.props.modelClass, model)
          }}>Delete</Button>
        </CardFooter>
      </React.Fragment>
    )
  }
}

export default ModelCard;