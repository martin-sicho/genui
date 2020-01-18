import React from "react";
import { ComponentWithObjects } from '../../../genui';
import ModelGrid from './ModelGrid';
// import { DropdownItem, DropdownMenu, DropdownToggle, UncontrolledDropdown } from 'reactstrap';

// function HeaderNav(props) {
//   return (<UncontrolledDropdown nav inNavbar>
//     <DropdownToggle nav caret>
//       Actions
//     </DropdownToggle>
//     <DropdownMenu right>
//       <UncontrolledDropdown>
//         <DropdownToggle nav>Add New...</DropdownToggle>
//         <DropdownMenu>
//           {
//             props.molSetChoices.map(choice =>
//               (<DropdownItem
//                 key={choice}
//                 onClick={() => {props.onMolSetChoice(choice, [])}}
//               >
//                 {choice}
//               </DropdownItem>)
//             )
//           }
//         </DropdownMenu>
//       </UncontrolledDropdown>
//     </DropdownMenu>
//   </UncontrolledDropdown>)
// }

class ModelsPage extends React.Component {

  componentDidMount() {
    // this.props.onHeaderChange(<HeaderNav {...this.props} molSetChoices={Object.keys(this.CLASS_TO_COMPONENT)} onMolSetChoice={this.props.handleAddMolSetList}/>);
  }

  render() {

    return (
      <div className="models-grid">
        <ComponentWithObjects
          {...this.props}
          objectListURL={new URL('models/', this.props.apiUrls.qsarRoot)}
          render={
            (models, handleAddModelList, handleAddModel, handleModelDelete) => {
              return <ModelGrid {...this.props} models={models} handleAddModel={handleAddModel} handleModelDelete={handleModelDelete} />
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