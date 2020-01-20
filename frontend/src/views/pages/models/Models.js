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
                onClick={() => {props.onModelAdd(choice, [])}}
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
      algorithmChoices : []
    }
  }

  handleAddNew = (model) => {
    this.setState({selectedToAdd : model})
  };

  fetchAlgorithms = () => {
    fetch(new URL('algorithms/', this.props.apiUrls.qsarRoot))
      .then(this.props.handleResponseErrors)
      .then((data) => {
        this.setState({ algorithmChoices: data })
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
  }

  render() {

    return (
      <div className="models-grid">
        <ComponentWithObjects
          {...this.props}
          objectListURL={new URL('models/', this.props.apiUrls.qsarRoot)}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid {...this.props} models={models} newModel={this.state.selectedToAdd} handleAddModel={handleAddModel} handleModelDelete={handleModelDelete} />
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