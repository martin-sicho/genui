import React from "react";
import { ComponentWithObjects } from '../../../genui';
import ModelGrid from './ModelGrid';
import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

function HeaderNav(props) {
  return (<UncontrolledDropdown nav inNavbar>
    <DropdownToggle nav caret>
      Actions
    </DropdownToggle>
    <DropdownMenu right>
      <UncontrolledDropdown>
        <DropdownToggle nav>Add New...</DropdownToggle>
        <DropdownMenu>
          {
            props.addChoices.map(choice =>
              (<DropdownItem
                key={choice.id}
                onClick={() => {props.onModelAdd(choice)}}
              >
                {choice.name}
              </DropdownItem>)
            )
          }
        </DropdownMenu>
      </UncontrolledDropdown>
    </DropdownMenu>
  </UncontrolledDropdown>)
}



class ModelsPage extends React.Component {

  constructor(props) {
    super(props);

    this.state = {
      selectedToAdd : null,
      algorithmChoices : [],
      descriptorChoices: [],
      metricsChoices: [],
      metricsLoaded : false,
      algorithmsLoaded : false,
      descriptorsLoaded: false,
    }
  }

  handleAddNew = (model) => {
    this.setState({selectedToAdd : model})
  };

  fetchAlgorithms = () => {
    fetch(new URL('algorithms/', this.props.apiUrls.qsarRoot))
      .then(this.props.handleResponseErrors)
      .then((data) => {
        this.setState({ algorithmChoices: data, algorithmsLoaded : true })
      })
    ;
  };

  fetchDescriptors = () => {
    fetch(new URL('descriptors/', this.props.apiUrls.qsarRoot))
      .then(this.props.handleResponseErrors)
      .then((data) => {
        this.setState({ descriptorChoices: data, descriptorsLoaded : true })
      })
    ;
  };

  fetchMetrics = () => {
    fetch(new URL('metrics/', this.props.apiUrls.qsarRoot))
      .then(this.props.handleResponseErrors)
      .then((data) => {
        this.setState({ metricsChoices: data, metricsLoaded : true })
      })
    ;
  };

  componentDidUpdate(prevProps, prevState, snapshot) {
    if (this.state.algorithmChoices && (prevState.algorithmChoices !== this.state.algorithmChoices)) {
      this.props.onHeaderChange(<HeaderNav {...this.props} addChoices={this.state.algorithmChoices} onModelAdd={this.handleAddNew}/>);
    }
  }

  componentDidMount() {
    this.fetchAlgorithms();
    this.fetchDescriptors();
    this.fetchMetrics();
  }

  render() {
    const {metricsLoaded, descriptorsLoaded, algorithmsLoaded} = {...this.state};
    if (!(metricsLoaded && descriptorsLoaded && algorithmsLoaded)) {
      return <div>Loading...</div>
    }

    return (
      <div className="models-grid">
        <ComponentWithObjects
          {...this.props}
          emptyClassName="QSARModel"
          objectListURL={new URL('models/', this.props.apiUrls.qsarRoot)}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid
                {...this.props}
                descriptors={this.state.descriptorChoices}
                metrics={this.state.metricsChoices}
                models={models}
                chosenAlgorithm={this.state.selectedToAdd}
                handleAddModel={
                  (...args) => {
                    this.setState({selectedToAdd : null});
                    return handleAddModel(...args)
                  }
                }
                handleModelDelete={handleModelDelete}
              />
            }
          }
        />
      </div>
    );
  }
}

function Models(props) {
  return (
    <ComponentWithObjects
      objectListURL={new URL('all/', props.apiUrls.compoundSetsRoot)}
      {...props}
      render={
        (
          compoundSets
        ) => {
          return (<ModelsPage
            {...props}
            compoundSets={compoundSets}
          />)
        }
      }
    />
  )
}

export default Models;