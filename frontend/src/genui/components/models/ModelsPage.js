import React from 'react';
import { ComponentWithObjects, ModelGrid } from '../../index';
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

    this.headerComponent = this.props.headerComponent ? this.props.headerComponent : HeaderNav;

    this.state = {
      selectedToAdd : null,
    }
  }

  handleAddNew = (model) => {
    this.setState({selectedToAdd : model})
  };

  componentDidMount() {
    const HeaderComp = this.headerComponent;
    this.props.onHeaderChange(<HeaderComp {...this.props} addChoices={this.props.algorithmChoices} onModelAdd={this.handleAddNew}/>);
  }

  render() {
    return (
      <div className="models-page">
        <ComponentWithObjects
          {...this.props}
          emptyClassName={this.props.modelClass}
          objectListURL={this.props.listURL}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid
                {...this.props}
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

export default ModelsPage;