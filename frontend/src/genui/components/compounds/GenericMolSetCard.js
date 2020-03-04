import React from "react";
import { Button, CardBody, CardFooter, CardHeader, CardSubtitle } from 'reactstrap';
import { TabWidget } from '../../index';

class GenericMolSetCard extends React.Component {
  abort = new AbortController();

  constructor(props) {
    super(props);
    this.created = new Date(this.props.molset.created);
    this.updated = new Date(this.props.molset.updated);
    this.molsetURL = new URL(`${this.props.molset.id}/`, this.props.molsetListUrl);
    this.moleculesURL = new URL(`${this.props.molset.id}/molecules/`, this.props.apiUrls.compoundSetsRoot);

    this.state = {
      molset : null
      , isUpdating : false
    }
  }

  handleDeleteSignal(deletedMolSet) {
    this.setState({isUpdating : true});
    if (this.props.hasOwnProperty("onMolsetDelete")) {
      this.props.onMolsetDelete(this.props.currentMolsetClass, deletedMolSet);
    }
  }

  componentDidMount() {
    fetch(this.molsetURL, {signal : this.abort.signal})
      .then(response => this.props.handleResponseErrors(response))
      .then(this.getMolSet)
      .catch(
        (error) => console.log(error)
      )
  }

  componentWillUnmount() {
    this.abort.abort();
  }

  getMolSet = (data) => {
    this.setState({
      molset : data
    })
  };

  updateMolSet = (data) => {
    const error_msg = 'Failed to update compound set from backend.';
    fetch(
      this.molsetURL
      , {
        method: 'PATCH'
        , body: JSON.stringify(data)
        , headers: {
          'Content-Type': 'application/json'
        }
      }
    )
      .then(response => this.props.handleResponseErrors(response, error_msg))
      .then(
        data => {
          this.setState({
            molset : data,
            isUpdating : true
          })
        }
      ).then(() => this.props.onTaskUpdate(
      () => this.setState({isUpdating : false})
      )
    )
      .catch(
        (error) => console.log(error)
      );
  };

  render() {
    const molset = this.state.molset;
    const isUpdating = this.state.isUpdating || this.props.tasksRunning;

    if (!molset) {
      return <div>Fetching...</div>
    }

    return (
      <React.Fragment>
        <CardHeader>{molset.name}</CardHeader>

        <CardBody className="scrollable">
          <CardSubtitle>
            <p>
              Created: {
              this.created.toLocaleDateString()
              + ' – ' + this.created.toLocaleTimeString()
            }
              <br/>
              Last Update: {
              this.updated.toLocaleDateString()
              + ' – ' + this.updated.toLocaleTimeString()
            }
            </p>
          </CardSubtitle>
          <TabWidget
            {...this.props}
            molsetIsUpdating={isUpdating}
            moleculesURL={this.moleculesURL}
            molset={molset}
          />
        </CardBody>

        <CardFooter>
          {/*TODO: transfer the update functionality to a separate Edit tab and use a form to do the update*/}
          <Button color="primary" disabled={isUpdating} onClick={() => this.updateMolSet({})}>{isUpdating ? 'Updating...' : 'Update Data'}</Button> <Button color="danger" disabled={isUpdating} onClick={() => {this.handleDeleteSignal(molset)}}>Delete</Button>
        </CardFooter>
      </React.Fragment>
    )
  }
}

export default GenericMolSetCard;